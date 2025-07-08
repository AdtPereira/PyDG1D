import numpy as np

from .dg3d_tools import *
from .dg1d_tools import jacobi_gauss
from .mesh3d import Mesh3D
from ..integrators.LSERK4 import *
from ..spatialDiscretization import *


class Maxwell3D(SpatialDiscretization):
    def __init__(self, n_order: int, mesh: Mesh3D, fluxType="Upwind"):
        assert n_order > 0
        assert mesh.number_of_elements() > 0

        self.n_order = n_order
        self.n_fp = self.number_of_nodes_per_face()
        self.n_faces = 4
        
        self.mesh = mesh
        self.fluxType = fluxType

        self.epsilon = np.ones(mesh.number_of_elements())
        self.mu = np.ones(mesh.number_of_elements())

        r, s, t = xyz_to_rst(*set_nodes_in_equilateral_tetrahedron(n_order))
        self.V = vandermonde(n_order, r, s, t)
        self.mass = massMatrix(n_order, r, s, t)
        self.Dr, self.Ds, self.Dt = derivateMatrix(n_order, r, s, t, self.V)
        self.x, self.y, self.z = nodesCoordinates(n_order, mesh)

        self.lift = lift(n_order)
        self.EToE, self.EToF = mesh.connectivityMatrices()

        self.rx, self.sx, self.tx, self.ry, self.sy, self.ty, self.rz, self.sz, self.tz, self.jacobian = geometricFactors(
            self.x, self.y, self.z, self.Dr, self.Ds, self.Dt)
        
        self.fmask, _, _, _, _ = buildFMask(n_order)

        self.nx, self.ny, self.nz, sJ = normals(self.x, self.y, self.z, self.Dr, self.Ds, self.Dt, n_order)
        self.f_scale = sJ/self.jacobian[self.fmask.ravel('F')]

        # Atributo para armazenar a fonte de excitação
        self.incident_wave_function = None
        self.current_corrections = {}

        # Dicionário para armazenar os mapas de nós nas interfaces entre grupos físicos
        self.buildMaps()

        # Build PML fields
        self.sgm_x = np.zeros_like(self.x)
        self.sgm_y = np.zeros_like(self.y)
        self.sgm_z = np.zeros_like(self.z)

    def set_current_source_corrections(self, corrections):
        """
        Armazena os campos da fonte calculados para o passo de tempo atual.
        """
        self.current_corrections = corrections

    def set_incident_wave_function(self, function_handle):
        """
        Armazena a função (method handle) que será usada para calcular a onda incidente.
        
        Parâmetros
        ----------
        function_handle : callable
            Uma função ou método que aceita (x, y, z, t) e retorna 6 componentes de campo.
        """
        self.incident_wave_function = function_handle

    def interpolate_dg_solution(self, uh, resolution=10):
        """
        Interpola a solução DG (uh) para uma grade de visualização densa dentro de cada elemento.

        Parâmetros
        ----------
        uh : ndarray
            Vetor da solução numérica com shape (Np, K), onde K é o número
            de elementos e Np o número de pontos por elemento.
        resolution : int, optional
            A ordem da grade de interpolação. Um valor maior gera uma visualização
            mais suave. O padrão é 10.

        Retorna
        -------
        x_interp, y_interp, z_interp : ndarray
            Arrays achatados com as coordenadas (x, y, z) dos pontos interpolados.
        uh_interp : ndarray
            Array achatado com os valores da solução interpolada nesses pontos.
        """
        # 1. Gerar uma grade densa de pontos no tetraedro de referência.
        # Usamos EquiNodes3D para criar pontos uniformemente espaçados em (r,s,t).
        r_interp, s_interp, t_interp = EquiNodes3D(resolution)

        # 2. Construir a matriz de interpolação.
        # `V` é a Vandermonde dos nós da solução. `V_interp` é a dos nós de interpolação.
        V_interp = vandermonde(self.n_order, r_interp, s_interp, t_interp) #
        invV = np.linalg.inv(self.V) # self.V foi calculado no __init__.
        interp_matrix = V_interp @ invV

        # 3. Mapear os pontos da grade de interpolação para as coordenadas físicas.
        # Usamos o método auxiliar definido acima.
        x_interp, y_interp, z_interp = map_reference_to_physical_nodes(self.mesh, r_interp, s_interp, t_interp)

        # 4. Aplicar a interpolação à solução `uh`.
        # O operador @ faz a multiplicação de matriz para cada elemento (coluna de uh).
        uh_interp = interp_matrix @ uh

        # 5. Retornar os arrays achatados, prontos para a plotagem com scatter.
        # A ordem 'F' (Fortran) é consistente com o resto do seu código.
        return (
            x_interp.ravel('F'),
            y_interp.ravel('F'),
            z_interp.ravel('F'),
            uh_interp.ravel('F'),
        )
    
    def compute_L2_error(self, uh, ua):
        """
        Calcula o erro global na norma L2 usando errᵗ diag(J) M err para cada elemento.

        Parâmetros
        ----------
        sp : DG1D
            Objeto com a discretização espacial.
        u_h : ndarray
            Solução numérica final do método DG (Np x K).
        ua : ndarray
            Solução analítica (Np x K).

        Retorno
        -------
        float
            Erro global na norma L2.
        """
        M = self.mass
        K = self.mesh.number_of_elements()
        err = ua - uh
        errL2 = np.zeros(K)

        for k in range(K):
            Jk = np.diag(self.jacobian[:, k])       # (Np x Np)
            ek = err[:, k][:, np.newaxis]           # (Np x 1)
            errL2[k] = (ek.T @ Jk @ M @ ek)[0, 0]   # (1 x 1)
        return np.sqrt(np.sum(errL2))

    def upml_sigma_fields(self, x, y, z, problem):
        """
        Avalia os campos σx(x) e σy(y) para as regiões de UPML.
        Parâmetros:
        - x, y: arrays de coordenadas onde avaliar os campos.
        - x_limit, y_limit: limites do domínio físico onde σ = 0.
        - sigma0_x, sigma0_y: coeficientes máximos de amortecimento.
        - p: expoente do perfil PML.
        Retorna:
        - sigma_x, sigma_y: arrays com os valores ponto a ponto para cada elemento.
        - dsigma_dx, dsigma_dy: derivadas dos campos σx e σy.
        """
        # --- Parâmetros Geométricos domínio ---
        Lx = problem['domain']['Lx']    # Dimensão total do domínio na direção x
        Ly = problem['domain']['Ly']    # Dimensão total do domínio na direção y
        Lz = problem['domain']['Lz']    # Dimensão total do domínio na direção z
        p = problem['pml']['pml_order'] # Ordem polinomial do perfil PML
        L = problem['pml']['L']         # Largura da camada da PML

        # --- Parâmetros PML ---
        x0 = Lx - L                     # Início do perfil PML na direção x
        y0 = Ly - L                     # Início do perfil PML na direção y
        z0 = Lz - L                     # Início do perfil PML na direção z

        # σx(x), σy(y) e σz(y)
        sigma_x = np.zeros_like(x)
        sigma_y = np.zeros_like(y)
        sigma_z = np.zeros_like(z)

        # Valor inicial e perfil PML
        sigma0_x = -np.log(problem['pml']['R'])
        sigma0_y = sigma0_x
        sigma0_z = sigma0_x        

        # Verifica se os limites são válidos
        if x0 <= 0 or y0 <= 0 or z0 <= 0:
            raise ValueError("Os limites x_limit, y_limit e z_limit devem ser positivos.")

        # σ^x(x)
        xPos = x >= x0
        xNeg = x <= -x0
        sigma_x[xPos] = sigma0_x * (x[xPos] - x0)**p
        sigma_x[xNeg] = sigma0_x * (x[xNeg] + x0)**p

        # σ^y(y)
        yPos = y >= y0
        yNeg = y <= -y0
        sigma_y[yPos] = sigma0_y * (y[yPos] - y0)**p
        sigma_y[yNeg] = sigma0_y * (y[yNeg] + y0)**p

        # σ^z(z)
        zPos = z >= z0
        zNeg = z <= -z0
        sigma_z[zPos] = sigma0_z * (z[zPos] - z0)**p
        sigma_z[zNeg] = sigma0_z * (z[zNeg] + z0)**p

        return sigma_x, sigma_y, sigma_z

    def buildMaps(self):
        '''
        Purpose: Build connectivity and boundary maps, and identify interface
                nodes and their indices between different physical groups.
        '''
        K = self.mesh.number_of_elements()
        Np = self.number_of_nodes_per_element()
        Nfaces = self.n_faces
        Nfp = self.n_fp

        # number volume nodes consecutively
        self.node_ids = np.reshape(np.arange(K*Np), [Np, K], 'F')
        vmapM = np.full([Nfp, Nfaces, K], 0)
        vmapP = np.full([Nfp, Nfaces, K], 0)
        # vmapI é agora um mapa de conectividade completo, com a mesma dimensão que os outros.
        vmapI = np.full([Nfp, Nfaces, K], 0, dtype=int)

        # Dicionário temporário para construir listas de nós de interface por pares de grupos
        tmp_interface_maps = {}

        # Encontra o índice dos nós da face em relação à ordem dos nós de volume
        for k1 in range(K):
            for f1 in range(Nfaces):
                vmapM[:, f1, k1] = self.node_ids[self.fmask[:, f1], k1]

        tmp = np.ones(Nfp)
        for k1 in range(K):
            for f1 in range(Nfaces):
                # Encontra o vizinho
                k2 = self.EToE[k1, f1]
                f2 = self.EToF[k1, f1]

                # Encontra os números dos nós de volume dos nós esquerdo e direito
                vidM = vmapM[:, f1, k1]
                vidP = vmapM[:, f2, k2]
                
                # Lógica de correspondência de nós para popular vmapP
                xM = np.outer(self.x.ravel('F')[vidM], tmp)
                yM = np.outer(self.y.ravel('F')[vidM], tmp)
                zM = np.outer(self.z.ravel('F')[vidM], tmp)
                
                xP = np.outer(self.x.ravel('F')[vidP], tmp)
                yP = np.outer(self.y.ravel('F')[vidP], tmp)
                zP = np.outer(self.z.ravel('F')[vidP], tmp)

                dist2 = (xM - xP.T)**2 + (yM - yP.T)**2 + (zM - zP.T)**2
                idM, idP = np.where(dist2 <= NODETOL)
                vmapP[idM, f1, k1] = vidP[idP]

                # --- NOVA LÓGICA PARA DEFINIR vmapI ---
                tagM = self.mesh.EToG[k1]
                tagP = self.mesh.EToG[k2]
                
                if tagM == tagP:
                    # Se os elementos estão na mesma região, vmapI é idêntico a vmapP.
                    # Isso cobre faces internas e também fronteiras de domínio (onde vmapP == vmapM).
                    vmapI[:, f1, k1] = vmapP[:, f1, k1]
                else: # Interface entre grupos físicos diferentes
                    # Se é uma interface física, vmapI é idêntico a vmapM.
                    vmapI[:, f1, k1] = vmapM[:, f1, k1]

        # Achata os mapas 3D para vetores 1D para uso eficiente
        self.vmapM = vmapM.ravel('F')
        self.vmapP = vmapP.ravel('F')
        self.vmapI = vmapI.ravel('F')

        # Identifica os nós de fronteira de domínio 
        # Uma fronteira é onde o vizinho é ele mesmo, logo vmapP == vmapM.
        boundary_mask = (self.vmapP == self.vmapM)
        self.mapB = np.where(boundary_mask)[0]
        self.vmapB = self.vmapM[boundary_mask]

        # 1. Criar uma máscara booleana para identificar as interfaces.
        # A máscara será 'True' onde vmapP e vmapI são diferentes.
        is_interface_face = self.vmapP != self.vmapI

        # vmapIM recebe os IDs dos nós do lado interno ('Minus') da interface.
        # Estes valores são extraídos de vmapM onde a máscara é verdadeira.
        self.vmapIM = self.vmapM[is_interface_face]

        # vmapIP recebe os IDs dos nós do lado externo ('Plus') da interface.
        # Estes valores são extraídos de vmapP onde a máscara é verdadeira.
        self.vmapIP = self.vmapP[is_interface_face]

        # 1. Obter os índices onde a condição de interface (vmapP != vmapI) é verdadeira.
        # A função np.where() retorna uma tupla de arrays; pegamos o primeiro elemento [0].
        self.mapI = np.where(is_interface_face)[0]
    
    def buildGroupInterfaceMaps(self, id_A, id_B):
        """
        Filtra os mapas de interface globais (vmapIM, vmapIP) para criar
        mapas de valores pareados para cada lado de uma interface específica.

        Para uma interface entre o Grupo A e o Grupo B, este método cria:
        - vmapIM_G{A}: Nós no lado 'Minus' que pertencem ao Grupo A.
        - vmapIP_G{A}: Nós parceiros no lado 'Plus' (pertencentes ao Grupo B).
        - vmapIM_G{B}: Nós no lado 'Minus' que pertencem ao Grupo B.
        - vmapIP_G{B}: Nós parceiros no lado 'Plus' (pertencentes ao Grupo A).

        Parâmetros
        ----------
        id_A : int
            O ID do primeiro grupo físico na interface.
        id_B : int
            O ID do segundo grupo físico na interface.
        """
        # 1. Determinar o elemento de origem ('Minus') para cada nó no mapa de interface.
        nodes_per_element_face = self.n_fp * self.n_faces
        element_indices_M = self.mapI // nodes_per_element_face

        # 2. Obter o grupo físico de cada elemento de origem.
        group_ids_M = self.mesh.EToG[element_indices_M]

        # 3. Criar máscaras para filtrar por grupo de origem.
        # Máscara para quando o lado 'Minus' da interface pertence ao Grupo A.
        mask_A_is_minus = (group_ids_M == id_A)
        # Máscara para quando o lado 'Minus' da interface pertence ao Grupo B.
        mask_B_is_minus = (group_ids_M == id_B)

        # 4. Aplicar as máscaras para criar os quatro mapas de valores desmembrados.
        
        # Par A: Visão da interface a partir do Grupo A
        vmapIM_G_A = self.vmapIM[mask_A_is_minus]
        vmapIP_G_A = self.vmapIP[mask_A_is_minus]

        # Par B: Visão da interface a partir do Grupo B
        vmapIM_G_B = self.vmapIM[mask_B_is_minus]
        vmapIP_G_B = self.vmapIP[mask_B_is_minus]

        # 5. Armazenar os novos mapas como atributos da classe usando setattr.
        setattr(self, f'vmapIM_G{id_A}', vmapIM_G_A)
        setattr(self, f'vmapIP_G{id_A}', vmapIP_G_A)        
        setattr(self, f'vmapIM_G{id_B}', vmapIM_G_B)
        setattr(self, f'vmapIP_G{id_B}', vmapIP_G_B)

        # --- NOVO Passo 6: Identificar e armazenar os elementos únicos de cada grupo ---
        elements_in_A = np.unique(element_indices_M[mask_A_is_minus])
        elements_in_B = np.unique(element_indices_M[mask_B_is_minus])        
        setattr(self, f'elements_G{id_A}', elements_in_A)
        setattr(self, f'elements_G{id_B}', elements_in_B)

        # Confirmação no console
        print(f"\nInterface value maps created for groups ({id_A}, {id_B}):")
        print(f" -> Found {len(elements_in_A)} unique elements for Group {id_A}.")
        print(f" -> Found {len(elements_in_B)} unique elements for Group {id_B}.")
        print(f" -> vmapIM_G{id_A} (Nodes in G{id_A}): {vmapIM_G_A.shape}")
        print(f" -> vmapIP_G{id_A} (Partners in G{id_B}): {vmapIP_G_A.shape}")
        print(f" -> vmapIM_G{id_B} (Nodes in G{id_B}): {vmapIM_G_B.shape}")
        print(f" -> vmapIP_G{id_B} (Partners in G{id_A}): {vmapIP_G_B.shape}")
        
        # Verificação de consistência:
        # Os nós do Grupo A devem ser os mesmos em vmapIM_G{A} e vmapIP_G{B}.
        # A ordenação pode ser diferente, então comparamos os conjuntos de nós.
        nodes_A_set1 = set(getattr(self, f'vmapIM_G{id_A}'))
        nodes_A_set2 = set(getattr(self, f'vmapIP_G{id_B}'))

        if nodes_A_set1 == nodes_A_set2:
            print(f" -> Node consistency check for Group {id_A} passed.")
        else:
            print(f" -> WARNING: Node consistency check for Group {id_A} failed.")

    def fieldsOnInterface(self, t):
        """
        Prepara os dicionários de correção para a interface TF/SF.

        Esta versão calcula os campos incidentes e os ordena de acordo com
        os mapas de cada lado da interface (SF->TF e TF->SF), retornando
        todos os dados necessários para a correção em computeJumps.
        """
        if not self.incident_wave_function:
            return {}

        # Obtenção dos mapas de índice para cada "lado" da interface
        try:
            map_sf = getattr(self, f'mapI_G{self.sf_group_id}')
            map_tf = getattr(self, f'mapI_G{self.tfz_group_id}')
        except AttributeError:
            # Se os mapas ainda não foram criados, retorna dicionário vazio
            return {}

        # --- Prepara correções para faces vistas do lado SF (SF -> TF) ---
        vmap_sf = self.vmapM[map_sf]
        x_sf = self.x.ravel('F')[vmap_sf]
        y_sf = self.y.ravel('F')[vmap_sf]
        z_sf = self.z.ravel('F')[vmap_sf]
        fields_inc_sf = self.incident_wave_function(x_sf, y_sf, z_sf, t)

        # --- Prepara correções para faces vistas do lado TF (TF -> SF) ---
        vmap_tf = self.vmapM[map_tf]
        x_tf = self.x.ravel('F')[vmap_tf]
        y_tf = self.y.ravel('F')[vmap_tf]
        z_tf = self.z.ravel('F')[vmap_tf]
        fields_inc_tf = self.incident_wave_function(x_tf, y_tf, z_tf, t)

        # Retorna uma estrutura de dados completa com todas as informações
        return {
            'map_scatter': map_sf,
            'map_total': map_tf,
            'fields_inc_scatter': fields_inc_sf, # (Hx_inc, Hy_inc, Hz_inc, Ex_inc, Ey_inc, Ez_inc)
            'fields_inc_total': fields_inc_tf
        }
    
    def get_minimum_node_distance(self):
        points, _ = jacobi_gauss(0, 0, self.n_order)
        return abs(points[0]-points[1])
    
    def get_dt_scale(self):
        """
        Calcula o fator de escala geométrico para o passo de tempo para cada 
        elemento tetraédrico 3D.

        Este fator é o raio da esfera inscrita no tetraedro (inradius), dado por
        r = 3V/A, onde V é o volume e A é a área de superfície total.
        Esta é uma métrica crucial para a condição CFL em 3D.
        """
        # Obtém os índices dos vértices para cada elemento (formato: [K, 4])
        # Onde K é o número de elementos.
        v = self.mesh.EToV

        # Obtém as coordenadas de todos os vértices da malha
        vx, vy, vz = self.mesh.vx, self.mesh.vy, self.mesh.vz

        # Reúne as coordenadas dos 4 vértices para cada um dos K elementos
        # p1, p2, p3, p4 terão formato [K, 3]
        p1 = np.vstack([vx[v[:, 0]], vy[v[:, 0]], vz[v[:, 0]]]).T
        p2 = np.vstack([vx[v[:, 1]], vy[v[:, 1]], vz[v[:, 1]]]).T
        p3 = np.vstack([vx[v[:, 2]], vy[v[:, 2]], vz[v[:, 2]]]).T
        p4 = np.vstack([vx[v[:, 3]], vy[v[:, 3]], vz[v[:, 3]]]).T

        # --- Calcula o Volume (V) ---
        # O volume de um tetraedro é V = |(p2-p1) . ((p3-p1) x (p4-p1))| / 6
        # Isso é equivalente a |det(M)| / 6, onde M é a matriz de vetores de aresta.
        mat = np.stack([p2 - p1, p3 - p1, p4 - p1], axis=1)
        volumes = np.abs(np.linalg.det(mat)) / 6.0

        # --- Calcula a Área de Superfície (A) ---
        # A área de um triângulo é 0.5 * ||(b-a) x (c-a)||

        # Face 1 (vértices p1, p2, p3)
        a1 = 0.5 * np.linalg.norm(np.cross(p2 - p1, p3 - p1), axis=1)
        # Face 2 (vértices p1, p2, p4)
        a2 = 0.5 * np.linalg.norm(np.cross(p2 - p1, p4 - p1), axis=1)
        # Face 3 (vértices p1, p3, p4)
        a3 = 0.5 * np.linalg.norm(np.cross(p3 - p1, p4 - p1), axis=1)
        # Face 4 (vértices p2, p3, p4)
        a4 = 0.5 * np.linalg.norm(np.cross(p3 - p2, p4 - p2), axis=1)

        surface_areas = a1 + a2 + a3 + a4

        # --- Calcula o Inradius (dt_scale) ---
        # Evita divisão por zero para elementos degenerados (área nula)
        dt_scale = np.zeros_like(surface_areas)
        non_zero_area = surface_areas > NODETOL
        dt_scale[non_zero_area] = 3.0 * volumes[non_zero_area] / surface_areas[non_zero_area]

        return dt_scale

    def get_mesh(self):
        return self.mesh

    def number_of_nodes_per_element(self):
        return int((self.n_order + 1) * (self.n_order + 2) * (self.n_order + 3) // 6)

    def number_of_nodes_per_face(self):
        return int((self.n_order + 1) * (self.n_order + 2) // 2)

    def count_boundary_elements(self):
        """
        Calcula o número de elementos que possuem pelo menos uma face na fronteira do domínio.

        Um elemento é considerado de fronteira se qualquer um de seus nós de face
        estiver contido no mapa de fronteira 'self.mapB'.

        Retorna
        -------
        int
            A quantidade total de elementos de fronteira únicos.
        """
        # Se mapB estiver vazio, não há nós de fronteira, logo, nenhum elemento de fronteira.
        if self.mapB.size == 0:
            return 0

        # O vetor vmapM é achatado com ordem 'F', onde o índice do elemento (k)
        # é o que varia mais lentamente. Podemos encontrá-lo com divisão inteira.
        # O divisor é o número total de pontos de face por elemento.
        nodes_per_element_face = self.n_fp * self.n_faces

        # Calcula o índice do elemento (k) para cada nó no mapa de fronteira.
        element_indices = self.mapB // nodes_per_element_face

        # Conta o número de elementos únicos na lista de índices resultante.
        num_boundary_elements = np.unique(element_indices).size

        return num_boundary_elements
    
    def buildEvolutionOperator(self):
        Np = self.number_of_nodes_per_element()
        K = self.mesh.number_of_elements()
        Nfields = 12            # Ex, Ey, Ez, Hx, Hy, Hz, Px, Py, Pz, Qx, Qy, Qz
        N = Nfields * Np * K    # total DOFs
        A = np.zeros((N,N))

        for i in range(N):
            fields = self.buildFields()
            node = i % Np
            elem = int(np.floor(i / Np)) % K
            field_block = i // (Np * K)

            # Set the appropriate field based on the block index
            if field_block == 0:        # Ex
                fields['Ex'][node, elem] = 1.0
            elif field_block == 1:      # Ey
                fields['Ey'][node, elem] = 1.0
            elif field_block == 2:      # Ez
                fields['Ez'][node, elem] = 1.0
            elif field_block == 3:      # Hx
                fields['Hx'][node, elem] = 1.0
            elif field_block == 4:      # Hy
                fields['Hy'][node, elem] = 1.0
            elif field_block == 5:      # Hz
                fields['Hz'][node, elem] = 1.0
            elif field_block == 6:      # Px
                fields['Px'][node, elem] = 1.0
            elif field_block == 7:      # Py
                fields['Py'][node, elem] = 1.0
            elif field_block == 8:      # Pz
                fields['Pz'][node, elem] = 1.0
            elif field_block == 9:      # Qx
                fields['Qx'][node, elem] = 1.0
            elif field_block == 10:     # Qy
                fields['Qy'][node, elem] = 1.0
            elif field_block == 11:     # Qz
                fields['Qz'][node, elem] = 1.0
            else:
                raise ValueError("Invalid field block index.")
            
            # Compute the right-hand side for the current fields
            fieldsRHS = self.computeRHS(fields)

            # Reshape the fields into a column vector and assign to the operator matrix
            q0 = np.vstack([
                fieldsRHS['Ex'].reshape(Np*K, 1, order='F'), 
                fieldsRHS['Ey'].reshape(Np*K, 1, order='F'), 
                fieldsRHS['Ez'].reshape(Np*K, 1, order='F'), 
                fieldsRHS['Hx'].reshape(Np*K, 1, order='F'),
                fieldsRHS['Hy'].reshape(Np*K, 1, order='F'),
                fieldsRHS['Hz'].reshape(Np*K, 1, order='F'),
                fieldsRHS['Px'].reshape(Np*K, 1, order='F'),
                fieldsRHS['Py'].reshape(Np*K, 1, order='F'),
                fieldsRHS['Pz'].reshape(Np*K, 1, order='F'),
                fieldsRHS['Qx'].reshape(Np*K, 1, order='F'),
                fieldsRHS['Qy'].reshape(Np*K, 1, order='F'),
                fieldsRHS['Qz'].reshape(Np*K, 1, order='F'),
            ])

            # Fill the operator matrix with the computed values
            A[:,i] = q0[:,0]

        return A

    def buildFields(self):
        # Electric fields components
        Ex = np.zeros([self.number_of_nodes_per_element(),
                       self.mesh.number_of_elements()])
        Ey = np.zeros(Ex.shape)
        Ez = np.zeros(Ex.shape)
        
        # Magnetic fields components
        Hx = np.zeros(Ex.shape)
        Hy = np.zeros(Ex.shape)
        Hz = np.zeros(Ex.shape)

        # UPML ADE E-fields
        Px = np.zeros(Ex.shape)
        Py = np.zeros(Ey.shape)
        Pz = np.zeros(Ez.shape)

        # UPML ADE H-fields
        Qx = np.zeros(Hx.shape)
        Qy = np.zeros(Hy.shape)
        Qz = np.zeros(Hz.shape)

        return {'Hx': Hx, 'Hy': Hy, 'Hz': Hz,
                'Ex': Ex, 'Ey': Ey, 'Ez': Ez,
                'Px': Px, 'Py': Py, 'Pz': Pz,
                'Qx': Qx, 'Qy': Qy, 'Qz': Qz}

    def computeZeroNormalFlux(self, dEx, dEy, dEz):
        f_Hx_zero = 0
        f_Hy_zero = 0
        f_Hz_zero = 0

        f_Ex_zero = 0
        f_Ey_zero = 0
        f_Ez_zero = 0

        if self.fluxType == "Upwind":
            ndotdE = self.nx * dEx + self.ny * dEy + self.nz * dEz
            f_Ex_zero -= dEx - ndotdE * self.nx
            f_Ey_zero -= dEy - ndotdE * self.ny
            f_Ez_zero -= dEz - ndotdE * self.nz

        elif self.fluxType == "Centered":
            pass

        else:
            raise ValueError("Invalid flux type.")

        return f_Hx_zero, f_Hy_zero, f_Hz_zero, f_Ex_zero, f_Ey_zero, f_Ez_zero

    def computeOneNormalFlux(self, dEx, dEy, dEz, dHx, dHy, dHz):

        f_Hx_one =  self.ny * dEz - self.nz * dEy
        f_Hy_one = -self.nx * dEz + self.nz * dEx
        f_Hz_one =  self.nx * dEy - self.ny * dEx

        f_Ex_one = -self.ny * dHz + self.nz * dHy
        f_Ey_one =  self.nx * dHz - self.nz * dHx
        f_Ez_one = -self.nx * dHy + self.ny * dHx

        return f_Hx_one, f_Hy_one, f_Hz_one, f_Ex_one, f_Ey_one, f_Ez_one

    def computeTwoNormalFlux(self, dHx, dHy, dHz):
        f_Hx_two = 0
        f_Hy_two = 0
        f_Hz_two = 0
        
        f_Ex_two = 0
        f_Ey_two = 0
        f_Ez_two = 0

        if self.fluxType == "Upwind":
            ndotdH = self.nx * dHx + self.ny * dHy + self.nz * dHz
            f_Hx_two += ndotdH * self.nx - dHx
            f_Hy_two += ndotdH * self.ny - dHy
            f_Hz_two += ndotdH * self.nz - dHz
        elif self.fluxType == "Centered":
            pass
        else:
            raise ValueError("Invalid flux type.")

        return f_Hx_two, f_Hy_two, f_Hz_two, f_Ex_two, f_Ey_two, f_Ez_two

    def computeFlux(self, Hx, Hy, Hz, Ex, Ey, Ez):
        dHx, dHy, dHz, dEx, dEy, dEz = self.computeJumps(Hx, Hy, Hz, Ex, Ey, Ez)
        
        f_Hx_zero, f_Hy_zero, f_Hz_zero, f_Ex_zero, f_Ey_zero, f_Ez_zero = self.computeZeroNormalFlux(dEx, dEy, dEz)
        f_Hx_one,   f_Hy_one,  f_Hz_one,  f_Ex_one,  f_Ey_one, f_Ez_one  = self.computeOneNormalFlux(dEx, dEy, dEz, dHx, dHy, dHz)
        f_Hx_two,   f_Hy_two,  f_Hz_two,  f_Ex_two, f_Ey_two,  f_Ez_two  = self.computeTwoNormalFlux(dHx, dHy, dHz)
        
        flux_Hx = f_Hx_zero + f_Hx_one + f_Hx_two
        flux_Hy = f_Hy_zero + f_Hy_one + f_Hy_two
        flux_Hz = f_Hz_zero + f_Hz_one + f_Hz_two

        flux_Ex = f_Ex_zero + f_Ex_one + f_Ex_two
        flux_Ey = f_Ey_zero + f_Ey_one + f_Ey_two
        flux_Ez = f_Ez_zero + f_Ez_one + f_Ez_two

        return flux_Hx, flux_Hy, flux_Hz, flux_Ex, flux_Ey, flux_Ez

    def fieldsOnBoundaryConditions(self, Hx, Hy, Hz, Ex, Ey, Ez):
        bcType = self.mesh.boundary_label
        if bcType == "PEC":
            Hbcx = Hx.transpose().take(self.vmapB)
            Hbcy = Hy.transpose().take(self.vmapB)
            Hbcz = Hz.transpose().take(self.vmapB)
            Ebcx = - Ex.transpose().take(self.vmapB)
            Ebcy = - Ey.transpose().take(self.vmapB)
            Ebcz = - Ez.transpose().take(self.vmapB)
        elif bcType == "PMC":
            Hbcx = - Hx.transpose().take(self.vmapB)
            Hbcy = - Hy.transpose().take(self.vmapB)
            Hbcz = - Hz.transpose().take(self.vmapB)
            Ebcx = Ex.transpose().take(self.vmapB)
            Ebcy = Ey.transpose().take(self.vmapB)
            Ebcz = Ez.transpose().take(self.vmapB)
        elif bcType == "SMA":
            Hbcx = Hx.transpose().take(self.vmapB) * 0.0
            Hbcy = Hx.transpose().take(self.vmapB) * 0.0
            Hbcz = Hz.transpose().take(self.vmapB) * 0.0
            Ebcx = Ex.transpose().take(self.vmapB) * 0.0
            Ebcy = Ey.transpose().take(self.vmapB) * 0.0
            Ebcz = Ez.transpose().take(self.vmapB) * 0.0
        elif bcType == "Periodic":
            Hbcx = Hx.transpose().take(self.vmapB[::-1])
            Hbcy = Hy.transpose().take(self.vmapB[::-1])
            Hbcz = Hz.transpose().take(self.vmapB[::-1])
            Ebcx = Ex.transpose().take(self.vmapB[::-1])
            Ebcy = Ey.transpose().take(self.vmapB[::-1])
            Ebcz = Ez.transpose().take(self.vmapB[::-1])
        else:
            raise ValueError("Invalid boundary label.")
        return Hbcx, Hbcy, Hbcz, Ebcx, Ebcy, Ebcz

    def computeJumps(self, Hx, Hy, Hz, Ex, Ey, Ez):
        Hbcx, Hbcy, Hbcz, Ebcx, Ebcy, Ebcz = self.fieldsOnBoundaryConditions(Hx, Hy, Hz, Ex, Ey, Ez)
        
        dHx = Hx.transpose().take(self.vmapM) - Hx.transpose().take(self.vmapP)
        dHy = Hy.transpose().take(self.vmapM) - Hy.transpose().take(self.vmapP)
        dHz = Hz.transpose().take(self.vmapM) - Hz.transpose().take(self.vmapP)
        
        dEx = Ex.transpose().take(self.vmapM) - Ex.transpose().take(self.vmapP)
        dEy = Ey.transpose().take(self.vmapM) - Ey.transpose().take(self.vmapP)
        dEz = Ez.transpose().take(self.vmapM) - Ez.transpose().take(self.vmapP)

        # Form field differences at faces
        dHx[self.mapB] = Hx.transpose().take(self.vmapB) - Hbcx
        dHy[self.mapB] = Hy.transpose().take(self.vmapB) - Hbcy
        dHz[self.mapB] = Hz.transpose().take(self.vmapB) - Hbcz
        
        dEx[self.mapB] = Ex.transpose().take(self.vmapB) - Ebcx
        dEy[self.mapB] = Ey.transpose().take(self.vmapB) - Ebcy
        dEz[self.mapB] = Ez.transpose().take(self.vmapB) - Ebcz

        # --- APLICAÇÃO DA CORREÇÃO TF/SF ---
        if self.current_corrections:
            # print("Aplicando correção TF/SF...")
            # Extrai os mapas e os campos incidentes da estrutura de dados
            map_scatter = self.current_corrections['map_scatter']
            map_total = self.current_corrections['map_total']
            
            # Campos incidentes ordenados para cada mapa
            # fields_inc é uma tupla (Hx, Hy, Hz, Ex, Ey, Ez)
            fields_inc_scatter = self.current_corrections['fields_inc_scatter']
            fields_inc_total = self.current_corrections['fields_inc_total']

            # Extrai os componentes de campo relevantes (Ey é o índice 4, Hz é o 2)
            Ey_inc_scatter = fields_inc_scatter[4]
            Hz_inc_scatter = fields_inc_scatter[2]            
            Ey_inc_total = fields_inc_total[4]
            Hz_inc_total = fields_inc_total[2]

            # 1. Aplica a correção para faces vistas do lado SF (SF(-) -> TF(+))
            # dEy = (Ey⁻ - Ey⁺) + Ey_inc⁺
            # Isso é equivalente a dEy = Ey⁻ - (Ey⁺ - Ey_inc⁺)
            dEy[map_scatter] = Ey.transpose().take(self.vmapM[map_scatter]) - (Ey.transpose().take(self.vmapP[map_scatter]) - Ey_inc_scatter)
            dHz[map_scatter] = Hz.transpose().take(self.vmapM[map_scatter]) - (Hz.transpose().take(self.vmapP[map_scatter]) - Hz_inc_scatter)

            # # 2. Aplica a correção para faces vistas do lado TF (TF(-) -> SF(+))
            # # dEy = (Ey⁻ - Ey⁺) - Ey_inc⁻
            # # Isso é equivalente a dEy = (Ey⁻ - Ey_inc⁻) - Ey⁺
            dEy[map_total] = Ey.transpose().take(self.vmapM[map_total]) - (Ey.transpose().take(self.vmapP[map_total]) + Ey_inc_total)
            dHz[map_total] = Hz.transpose().take(self.vmapM[map_total]) - (Hz.transpose().take(self.vmapP[map_total]) + Hz_inc_total)

        dHx = dHx.reshape(self.n_fp*self.n_faces, self.mesh.number_of_elements(), order='F')
        dHy = dHy.reshape(self.n_fp*self.n_faces, self.mesh.number_of_elements(), order='F')
        dHz = dHz.reshape(self.n_fp*self.n_faces, self.mesh.number_of_elements(), order='F')
        
        dEx = dEx.reshape(self.n_fp*self.n_faces, self.mesh.number_of_elements(), order='F')
        dEy = dEy.reshape(self.n_fp*self.n_faces, self.mesh.number_of_elements(), order='F')
        dEz = dEz.reshape(self.n_fp*self.n_faces, self.mesh.number_of_elements(), order='F')

        return dHx, dHy, dHz, dEx, dEy, dEz

    def computeRHS(self, fields):
        # Compute the right-hand side of Maxwell's equations
        Hx, Hy, Hz = fields['Hx'], fields['Hy'], fields['Hz']
        Ex, Ey, Ez = fields['Ex'], fields['Ey'], fields['Ez']

        # Compute the right-hand side of ADE equations
        Px, Py, Pz = fields['Px'], fields['Py'], fields['Pz']
        Qx, Qy, Qz = fields['Qx'], fields['Qy'], fields['Qz']

        # Compute the RHS for the time evolution
        flux_Hx, flux_Hy, flux_Hz, flux_Ex, flux_Ey, flux_Ez = self.computeFlux(
            Hx, Hy, Hz, Ex, Ey, Ez)

        # evaluate local spatial derivatives
        curlHx, curlHy, curlHz = Curl3D(
            Hx, Hy, Hz,
            self.Dr, self.Ds, self.Dt,
            self.rx, self.sx, self.tx,
            self.ry, self.sy, self.ty,
            self.rz, self.sz, self.tz
        )

        curlEx, curlEy, curlEz = Curl3D(
            Ex, Ey, Ez,
            self.Dr, self.Ds, self.Dt,
            self.rx, self.sx, self.tx,
            self.ry, self.sy, self.ty,
            self.rz, self.sz, self.tz
        )       
        
        # calculate Maxwell's right hand side
        # missing material epsilon/mu
        rhs_Hx = -curlEx + np.matmul(self.lift, self.f_scale * flux_Hx)/2.0
        rhs_Hy = -curlEy + np.matmul(self.lift, self.f_scale * flux_Hy)/2.0
        rhs_Hz = -curlEz + np.matmul(self.lift, self.f_scale * flux_Hz)/2.0

        rhs_Ex = +curlHx + np.matmul(self.lift, self.f_scale * flux_Ex)/2.0
        rhs_Ey = +curlHy + np.matmul(self.lift, self.f_scale * flux_Ey)/2.0
        rhs_Ez = +curlHz + np.matmul(self.lift, self.f_scale * flux_Ez)/2.0

        # -------- Termos extras da formulação PML  ----------
        # -------- missing material epsilon/mu      ----------

        # UPML E-fields
        rhs_Ex += - Px + (+self.sgm_x - self.sgm_y - self.sgm_z) * Ex
        rhs_Ey += - Py + (-self.sgm_x + self.sgm_y - self.sgm_z) * Ey    
        rhs_Ez += - Pz + (-self.sgm_x - self.sgm_y + self.sgm_z) * Ez

        # UPML H-fields
        rhs_Hx += - Qx + (+self.sgm_x - self.sgm_y - self.sgm_z) * Hx
        rhs_Hy += - Qy + (-self.sgm_x + self.sgm_y - self.sgm_z) * Hy
        rhs_Hz += - Qz + (-self.sgm_x - self.sgm_y + self.sgm_z) * Hz

        # UPML ADE fields
        sgm_xx = self.sgm_x ** 2
        sgm_yy = self.sgm_y ** 2
        sgm_zz = self.sgm_z ** 2

        sgm_xy = self.sgm_x * self.sgm_y
        sgm_xz = self.sgm_x * self.sgm_z
        sgm_yz = self.sgm_y * self.sgm_z

        rhs_Px = (sgm_xx - sgm_xy - sgm_xz + sgm_yz) * Ex - self.sgm_x * Px
        rhs_Py = (sgm_yy - sgm_xy + sgm_xz - sgm_yz) * Ey - self.sgm_y * Py
        rhs_Pz = (sgm_zz + sgm_xy - sgm_xz - sgm_yz) * Ez - self.sgm_z * Pz

        rhs_Qx = (sgm_xx - sgm_xy - sgm_xz + sgm_yz) * Hx - self.sgm_x * Qx
        rhs_Qy = (sgm_yy - sgm_xy + sgm_xz - sgm_yz) * Hy - self.sgm_y * Qy
        rhs_Qz = (sgm_zz + sgm_xy - sgm_xz - sgm_yz) * Hz - self.sgm_z * Qz

        return {'Hx': rhs_Hx, 'Hy': rhs_Hy, 'Hz': rhs_Hz,
                'Ex': rhs_Ex, 'Ey': rhs_Ey, 'Ez': rhs_Ez,
                'Px': rhs_Px, 'Py': rhs_Py, 'Pz': rhs_Pz,
                'Qx': rhs_Qx, 'Qy': rhs_Qy, 'Qz': rhs_Qz}
