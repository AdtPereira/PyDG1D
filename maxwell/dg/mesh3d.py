import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

N_FACES = 4

class Mesh3D:
    def __init__(self, vx, vy, vz, EToV, boundary_label="PEC"):
        assert vx.shape == vy.shape == vz.shape
        assert np.max(np.max(EToV))+1  == vx.shape[0]
        
        self.dimension = 3
        self.vx = vx
        self.vy = vy
        self.vz = vz
        self.EToV = EToV        
        self.boundary_label = boundary_label

        # Vetor que mapeia cada elemento (índice) a uma tag de grupo (valor)
        self.EToG = np.zeros(self.number_of_elements(), dtype=int)

        # Dicionário auxiliar para armazenar nomes dos grupos
        self.group_names = {0: 'Default'}


    def number_of_vertices(self):
        return self.vx.shape[0]


    def number_of_elements(self):
        return self.EToV.shape[0]
    
    
    def connectivityMatrices(self):
        '''
            % function [EToE,EToF]= tiConnect3D(EToV)
            % Purpose: tetrahedral face connect algorithm due to Toby Isaac

            Conecta faces de elementos tetraédricos.

            Parâmetros:
            -----------
            EToV : ndarray (K x 4)
                Matriz de conectividade dos elementos tetraédricos, onde cada linha contém os índices dos 4 vértices do elemento.

            Retorna:
            --------
            EToE : ndarray (K x 4)
                Matriz de conectividade elemento a elemento.
            EToF : ndarray (K x 4)
                Matriz de conectividade elemento a face.
        '''
        K = self.number_of_elements()
        Nnodes = self.number_of_vertices()
                
        # create list of all faces 1, then 2, 3 & 4
        fnodes = np.concatenate(
            (self.EToV[:,[0, 1, 2]],    # face oposta ao vértice 3
             self.EToV[:,[0, 1, 3]],    # face oposta ao vértice 2
             self.EToV[:,[1, 2, 3]],    # face oposta ao vértice 0  
             self.EToV[:,[0, 2, 3]])    # face oposta ao vértice 1
        )

        # ordena vértices de cada face para identificação única
        fnodes = np.sort(fnodes, axis=1)

        # set up default element to element and Element to faces connectivity
        EToE = np.outer(np.arange(K, dtype=int), np.ones((N_FACES,), dtype=int))
        EToF = np.outer(np.ones((K,), dtype=int), np.arange(N_FACES, dtype=int))

        # Gera identificadores únicos por face
        ids = fnodes[:, 0] * Nnodes * Nnodes + fnodes[:, 1] * Nnodes + fnodes[:, 2]

        spNodeToNode = np.concatenate([
            ids.reshape(-1, 1),
            np.arange(N_FACES * K).reshape(-1, 1),
            EToE.reshape(-1, 1, order='F'),
            EToF.reshape(-1, 1, order='F')
        ], axis=1)

        spNodeToNode = spNodeToNode[np.argsort(spNodeToNode[:, 0])]

        # Encontrar faces coincidentes
        matches = np.where(spNodeToNode[:-1, 0] == spNodeToNode[1:, 0])[0]
        matchL = np.concatenate([spNodeToNode[matches],     spNodeToNode[matches + 1]])
        matchR = np.concatenate([spNodeToNode[matches + 1], spNodeToNode[matches]])

        # Atualiza conectividades
        EToE.T.flat[matchL[:, 1]] = matchR[:, 2]
        EToF.T.flat[matchL[:, 1]] = matchR[:, 3]

        return EToE, EToF


    def get_elements_by_group(self, group_tag):
        """
        Retorna os índices de todos os elementos que pertencem a um dado grupo.
        """
        # np.where retorna uma tupla; pegamos o primeiro elemento
        return np.where(self.EToG == group_tag)[0]


    def box_physical_group(self, group_tag, group_name, half_dim, center=(0, 0, 0), group_info=True):
        """
        Atribui uma nova tag de grupo aos elementos contidos em uma caixa.
        Modifica o vetor self.EtoG diretamente.
        """
        if group_tag == 0:
            print("Erro: A tag 0 é reservada para o grupo 'Default'. Use outra tag.")
            return np.array([], dtype=int)
            
        print(f"Atribuindo grupo '{group_name}' (Tag: {group_tag})...")
        min_coords = np.array(center) - half_dim
        max_coords = np.array(center) + half_dim
        
        elements_to_assign = []
        for k in range(self.number_of_elements()):
            vert_indices = self.EToV[k, :]
            v_coords_x = self.vx[vert_indices]
            v_coords_y = self.vy[vert_indices]
            v_coords_z = self.vz[vert_indices]

            all_vertices_in_box = (
                np.all((v_coords_x >= min_coords[0]) & (v_coords_x <= max_coords[0])) and
                np.all((v_coords_y >= min_coords[1]) & (v_coords_y <= max_coords[1])) and
                np.all((v_coords_z >= min_coords[2]) & (v_coords_z <= max_coords[2]))
            )

            if all_vertices_in_box:
                elements_to_assign.append(k)

        # Converte para array para indexação avançada
        elements_to_assign = np.array(elements_to_assign, dtype=int)
        
        if elements_to_assign.size > 0:
            # Atribui a nova tag a todos os elementos encontrados de uma vez
            self.EToG[elements_to_assign] = group_tag
        
        # Armazena o nome do novo grupo
        self.group_names[group_tag] = group_name
        
        print(f"-> {len(elements_to_assign)} elementos foram atribuídos à tag {group_tag}.")

        if group_info:
            # 3. Relatório final usando os novos métodos de acesso
            print("\n--- Relatório Final de Grupos Físicos ---")
            print("Grupos definidos:", self.group_names)

            for tag, name in self.group_names.items():
                # Usa o novo método para buscar os elementos por grupo
                elements_in_group = self.get_elements_by_group(tag)
                print(f"  - Grupo '{name}' (Tag {tag}): {len(elements_in_group)} elementos.")

            assigned_tets = len(self.get_elements_by_group(group_tag))
            unassigned_tets_count = len(self.get_elements_by_group(0))

            assert self.number_of_elements() == assigned_tets + unassigned_tets_count, \
            f"Erro: Total de elementos ({self.number_of_elements()}) não corresponde à soma dos atribuídos ({assigned_tets}) e não atribuídos ({unassigned_tets_count})."
        
        return elements_to_assign
    

    def plot_mesh(self, title="Malha 3D", show_vertices=False, alpha=0.25):
        """
        Plota a malha 3D de tetraedros usando matplotlib.

        Parâmetros:
        -----------
        title : str
            O título do gráfico.
        show_vertices : bool
            Se True, plota os vértices como pontos.
        alpha : float
            Nível de transparência das faces do tetraedro (entre 0 e 1).
        """
        if self.number_of_elements() == 0:
            print("A malha não contém elementos para plotar.")
            return

        fig = plt.figure(figsize=(10, 8))
        ax = fig.add_subplot(111, projection='3d')

        # Coleta todas as faces de todos os tetraedros
        all_faces = []
        for k in range(self.number_of_elements()):
            # Índices dos 4 vértices do tetraedro atual
            vert_indices = self.EToV[k, :]
            
            # Coordenadas dos 4 vértices
            v = np.array([
                self.vx[vert_indices],
                self.vy[vert_indices],
                self.vz[vert_indices]
            ]).T

            # Define as 4 faces triangulares do tetraedro
            # Face 1: vértices 0, 1, 2
            # Face 2: vértices 0, 1, 3
            # Face 3: vértices 0, 2, 3
            # Face 4: vértices 1, 2, 3
            all_faces.append([v[0], v[1], v[2]])
            all_faces.append([v[0], v[1], v[3]])
            all_faces.append([v[0], v[2], v[3]])
            all_faces.append([v[1], v[2], v[3]])

        # Cria a coleção de polígonos para uma renderização eficiente
        collection = Poly3DCollection(
            all_faces, 
            alpha=alpha, 
            facecolors='cyan', 
            linewidths=1, 
            edgecolors='k'
        )
        ax.add_collection3d(collection)

        # Plota os vértices se solicitado
        if show_vertices:
            ax.scatter(self.vx, self.vy, self.vz, color='r', s=10, label='Vértices')
            ax.legend()
            
        # Ajusta os limites do gráfico para garantir que a malha inteira seja visível
        ax.set_xlim(self.vx.min(), self.vx.max())
        ax.set_ylim(self.vy.min(), self.vy.max())
        ax.set_zlim(self.vz.min(), self.vz.max())
        ax.set_box_aspect([
            np.ptp(self.vx), 
            np.ptp(self.vy), 
            np.ptp(self.vz)
        ])  # Proporção de aspecto para evitar distorções

        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.set_zlabel('Z')
        ax.set_title(title)
        
        plt.show()
    

