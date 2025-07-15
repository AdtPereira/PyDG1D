import numpy as np
import gmsh
import matplotlib.tri as mtri

N_FACES = 3

class Mesh2D:
    def __init__(self, vx, vy, EToV, physical_groups = None, boundary_label="PEC"):
        assert vx.shape == vy.shape
        assert np.max(np.max(EToV)) + 1  == vx.shape[0]
        
        self.dimension = 2
        self.vx = vx
        self.vy = vy
        self.EToV = EToV        
        self.boundary_label = boundary_label        

        # Dicionário para armazenar grupos físicos
        self.physical_groups = physical_groups if physical_groups is not None else {}

        # Gera e armazena os mapeamentos de grupos físicos
        self.EToG = self.connectivityElementsToPhysicalGroups()
        self.FToG = self.connectivityFacesToPhysicalGroups()


    def number_of_vertices(self):
        return self.vx.shape[0]


    def number_of_elements(self):
        return self.EToV.shape[0]


    def get_elements_by_group(self, group_tag: int) -> np.ndarray:
        """ 
        Retorna os índices de todos os elementos que pertencem a um dado grupo.
        
        Args:
            group_tag (int): A tag do grupo físico a ser buscado.

        Returns:
            np.ndarray: Um array com os índices dos elementos que pertencem ao grupo.
        """
        # np.where encontra os índices onde a condição (EToG == group_tag) é verdadeira.
        # Ele retorna uma tupla de arrays, então pegamos o primeiro elemento [0].
        return np.where(self.EToG == group_tag)[0]
    
    
    def getTriangulation(self):
        return mtri.Triangulation(self.vx, self.vy, self.EToV.tolist())
    

    def connectivityMatrices(self):
        '''
            function [EToE,EToF]= connectivity(self)
            Purpose: triangle face connect algorithm due to Toby Isaac
        '''
        K = self.number_of_elements()
        Nnodes = self.number_of_vertices()
                
        # create list of all faces 0, then 1, & 2
        fnodes = np.concatenate(
            (self.EToV[:,[0,1]],    # Todas as faces 0 (vértices 0 e 1 de cada elemento)
             self.EToV[:,[1,2]],    # Todas as faces 1 (vértices 1 e 2 de cada elemento)    
             self.EToV[:,[2,0]])    # Todas as faces 2 (vértices 2 e 0 de cada elemento)
        )
        fnodes = np.sort(fnodes, 1)

        # set up default element to element and Element to faces connectivity
        EToE = np.outer(np.arange(0, K, 1, dtype=int), np.ones((1,N_FACES), dtype=int))
        EToF = np.outer(np.ones((K,1), dtype=int),     np.arange(0, N_FACES, 1, dtype=int))

        # uniquely number each set of three faces by their node numbers 
        id = fnodes[:,0] * Nnodes + fnodes[:,1] + 1
        spNodeToNode = np.concatenate((
            id.reshape(id.size, 1), 
            np.arange(0, N_FACES*K, 1, dtype=int).reshape(N_FACES*K, 1), 
            EToE.reshape(N_FACES*K, 1, order='F'), 
            EToF.reshape(N_FACES*K, 1, order='F')
        ), 1)
        spNodeToNode = np.int64(spNodeToNode)

        # Now we sort by global face number.
        sorted_nodes = spNodeToNode[np.argsort(spNodeToNode[:,0]), :]

        # find matches in the sorted face list
        indices = np.where(sorted_nodes[:-1,0]==sorted_nodes[1:, 0])[0]

        # make links reflexive 
        matchL = np.concatenate((sorted_nodes[indices,  :], sorted_nodes[indices+1,:]))
        matchR = np.concatenate((sorted_nodes[indices+1,:], sorted_nodes[indices,  :]))

        # insert matches
        EToE.transpose().flat[matchL[:,1]] = matchR[:,2]
        EToF.transpose().flat[matchL[:,1]] = matchR[:,3]

        return EToE, EToF


    def connectivityElementsToPhysicalGroups(self):
        """
        Constrói o mapeamento EToG (Elemento para Grupo Físico).

        Cada posição 'k' no array EToG armazena a tag do grupo físico
        à qual o elemento 'k' da malha global (self.EToV[k]) pertence.

        A lógica se baseia na criação de um mapa de busca onde a chave é uma
        representação canônica de um elemento (seus vértices ordenados) e o valor
        é seu índice global 'k'. O método então percorre os grupos físicos
        de dimensão 2 e preenche o array EToG.

        Retorna:
            np.ndarray: Um vetor de inteiros de tamanho K (número de elementos),
                        onde EToG[k] é a tag do grupo físico do elemento k.
        """
        K = self.number_of_elements()
        # Inicializa EToG com 0 ou outra tag padrão.
        EToG = np.zeros(K, dtype=int)

        if not self.physical_groups:
            return EToG

        # 1. Cria um mapa de busca: (tupla de vértices ordenados) -> índice do elemento global (k)
        # Uma tupla de vértices ordenados é um identificador único e canônico para um elemento.
        element_map = {tuple(sorted(self.EToV[k])): k for k in range(K)}

        # 2. Itera sobre todos os grupos físicos fornecidos
        for tag, group_data in self.physical_groups.items():
            # Considera apenas grupos de elementos de superfície (dimensão 2)
            if group_data.get('dim') == self.dimension:
                group_EToV = group_data.get('EToV', [])
                
                # 3. Para cada elemento no grupo, encontra seu índice global e atribui a tag
                for element_vertices in group_EToV:
                    # Cria a chave canônica para o elemento do grupo
                    key = tuple(sorted(element_vertices))
                    
                    # Procura o elemento no mapa da malha global
                    if key in element_map:
                        # Atribui a tag do grupo físico à posição correta em EToG
                        EToG[element_map[key]] = tag
        
        return EToG


    def connectivityFacesToPhysicalGroups(self):
        """
        Constrói o mapeamento FToG (Face para Grupo Físico) com uma estrutura
        matricial (K, N_FACES), análoga a EToE e EToF.
        """
        K = self.number_of_elements()
        FToG = np.zeros((K, N_FACES), dtype=int)

        if not self.physical_groups:
            return FToG

        # 1. Mapeia cada face geométrica para sua(s) coordenada(s) (k, f)
        face_map = {}
        for k in range(K):
            vertex = self.EToV[k]
            local_faces = [[vertex[0], vertex[1]],
                           [vertex[1], vertex[2]],
                           [vertex[2], vertex[0]]]

            for f, face_vertex in enumerate(local_faces):
                key = tuple(sorted(face_vertex))
                if key not in face_map:
                    face_map[key] = []
                face_map[key].append((k, f))

        # 2. Itera sobre os grupos físicos de dimensão 1 (curvas/contornos)
        for tag, group_data in self.physical_groups.items():
            if group_data.get('dim') == self.dimension - 1:
                boundary_faces = group_data.get('EToV', [])
                
                for face_vertex in boundary_faces:
                    key = tuple(sorted(face_vertex))
                    
                    if key in face_map:
                        coords_list = face_map[key]
                        
                        # --- CORREÇÃO APLICADA AQUI ---
                        # Itera sobre todos os pares (k, f) que compartilham esta face.
                        # Para uma face externa, a lista terá 1 item.
                        # Para uma face interna, a lista terá 2 itens.
                        for k, f in coords_list:
                            FToG[k, f] = tag
                            
        return FToG


class GmshMeshReader:
    """
    Uma classe para ler arquivos de malha .msh gerados pelo Gmsh.

    Esta classe funciona como um gerenciador de contexto para garantir
    que os recursos do Gmsh sejam inicializados e finalizados corretamente.

    Processa a malha para extrair dados (`vx`, `vy`, `EToV`) para cada
    grupo físico definido no arquivo.
    """
    def __init__(self, filename: str, problem_dim: int):
        self.filename = filename
        self.ProblemDim = problem_dim

        # Atributos para armazenar os dados da malha global
        self.vx = None
        self.vy = None
        self.EToV = None
        # self.physical_groups = None


    def __enter__(self):
        gmsh.initialize()
        gmsh.open(self.filename)
        return self


    def __exit__(self, exc_type, exc_val, exc_tb):
        gmsh.finalize()


    def read_global_mesh(self) -> None:
        """
        Lê a malha 2D global completa do arquivo, combinando todos os elementos
        dos grupos físicos de dimensão 2. Este método reutiliza os métodos
        auxiliares para construir a malha global.
        """
        # 1. Identificar todos os grupos físicos de dimensão 2
        all_physical_groups = gmsh.model.getPhysicalGroups()
        dim2_groups = [(dim, tag) for dim, tag in all_physical_groups if dim == self.ProblemDim]

        if not dim2_groups:
            raise RuntimeError(f"Nenhum grupo físico de dimensão {self.ProblemDim} foi encontrado na malha.")

        # 2. Coletar dados de nós de todos os grupos 2D para criar um conjunto único
        all_nodes_data = {}
        for dim, tag in dim2_groups:
            # Chama o método auxiliar para cada grupo 2D
            group_info = self._getPhysicalNodesDict((dim, tag))
            # Adiciona/atualiza o dicionário global com os nós de cada grupo
            all_nodes_data.update(group_info['nodes'])

        # 3. Extrair vértices únicos e criar o mapa de 'tag' para 'índice'
        # A função _extractVertices é efetivamente replicada aqui para o conjunto global
        sorted_nodes = sorted(all_nodes_data.items())
        node_tags = [key for key, _ in sorted_nodes]
        self.vx = np.array([value[0] for _, value in sorted_nodes])
        self.vy = np.array([value[1] for _, value in sorted_nodes])
        
        # Este mapa é a chave para a re-indexação
        tag_to_idx = {tag: i for i, tag in enumerate(node_tags)}

        # 4. Coletar a conectividade (EToV) de cada grupo 2D usando _getEToV
        list_of_EToV_with_tags = []
        for dim, tag in dim2_groups:
            # Chamamos _getEToV com index_based=1 para que ele retorne as tags
            # originais do GMSH, em vez de tentar convertê-las para base 0.
            EToV_group = self._getEToV((dim, tag), index_based=1)
            list_of_EToV_with_tags.append(EToV_group)

        # 5. Concatenar todas as matrizes EToV
        if not list_of_EToV_with_tags:
            raise RuntimeError("Não foi possível extrair a conectividade de nenhum grupo 2D.")
        
        EToV_with_tags = np.concatenate(list_of_EToV_with_tags, axis=0)

        # 6. Re-indexar a matriz EToV global usando o mapa tag_to_idx
        # np.vectorize aplica a busca no dicionário a cada elemento de forma eficiente
        self.EToV = np.vectorize(tag_to_idx.get)(EToV_with_tags).astype(int)


    def view(self):
        """
        Abre a interface gráfica do Gmsh para visualizar a malha.
        """
        gmsh.fltk.run()


    def _getEToV(self, GroupDimTag: tuple[int, int], index_based=0) -> np.ndarray:
        
        EToV = []
        GroupEntitiesTags = gmsh.model.getEntitiesForPhysicalGroup(*GroupDimTag)

        for EntityTag in GroupEntitiesTags:
            elemTypes, _, elemNodeTags = gmsh.model.mesh.getElements(GroupDimTag[0], EntityTag)

            for elemType, elemNodes in zip(elemTypes, elemNodeTags):
                _, _, _, numNodes, _, _ = gmsh.model.mesh.getElementProperties(elemType)
                
                for i in range(len(elemNodes) // numNodes):
                    conn = elemNodes[numNodes * i: numNodes * (i + 1)]
                    EToV.append([int(tag) for tag in conn])
                
        # Se index_based for False, converter os índices para base 0
        if not index_based:
            EToV = [np.array(conn) - 1 for conn in EToV]

        return np.array(EToV, dtype=int)


    def _getPhysicalNodesDict(self, GroupDimTag: tuple[int, int]) -> dict:
        """
        Retorna os dados de nós para um único grupo físico especificado.

        Args:
            GroupDimTag (tuple[int, int]): A tupla (dimensão, tag) do grupo físico.

        Returns:
            dict: Um dicionário contendo o nome e os dados dos nós do grupo.
        """
        dim, GroupTag = GroupDimTag
        NodeTags, NodeCoords = gmsh.model.mesh.getNodesForPhysicalGroup(dim, GroupTag)
        NodeTags = [int(tag) for tag in NodeTags.astype(np.uint64)]

        group_data = {
            "name": gmsh.model.getPhysicalName(dim, GroupTag),
            "nodes": {
                int(node): (
                    float(NodeCoords[3*i]),
                    float(NodeCoords[3*i + 1]),
                    float(NodeCoords[3*i + 2]),
                ) for i, node in enumerate(NodeTags)
            }
        }
        return group_data
    

    def _extractVertices(self, nodesData):
        sorted_nodes = sorted(nodesData.items())
        vx = np.array([value[0] for key_, value in sorted_nodes])
        vy = np.array([value[1] for key_, value in sorted_nodes])
        return vx, vy


    def get_physical_group_data(self) -> dict:
        """
        Coleta e estrutura os dados de malha para TODOS os grupos físicos (1D, 2D, etc.).

        Retorna:
            dict: Um dicionário onde as chaves são os nomes dos grupos físicos
                  e os valores são outros dicionários contendo {'dim', 'vx', 'vy', 'EToV'}
                  para cada grupo.

        Exemplo de retorno:
        {
            "outerBoundary": {
                "dim": 1,
                "vx": [x0, x1, ...],
                "vy": [y0, y1, ...],
                "EToV": [[0, 1], [1, 2], ...]
            },
            "innerDomain": {
                "dim": 2,
                "vx": [x0, x1, ...],
                "vy": [y0, y1, ...],
                "EToV": [[0, 1, 5], [1, 2, 6], ...]
            }
        }
        """
        physical_group_dict = {}

        # Itera sobre a lista de grupos físicos, que já contém a dimensão e a tag
        for dimGroup, tagGroup in gmsh.model.getPhysicalGroups():
            # Chama o método padronizado para obter os dados do grupo atual
            group_info = self._getPhysicalNodesDict((dimGroup, tagGroup))
            group_name = group_info['name']
            nodes_data = group_info['nodes']

            # Extrai os vértices e a conectividade para este grupo específico
            # A função _getEToV funciona para qualquer dimensão
            EToV = self._getEToV((dimGroup, tagGroup))
            vx, vy = self._extractVertices(nodes_data)

            # Armazena os dados, incluindo a dimensão do grupo
            physical_group_dict[tagGroup] = {
                'name': group_name,
                'dim': dimGroup,
                'vx': vx,
                'vy': vy,
                'EToV': EToV
            }
        
        return physical_group_dict
    

def readFromGambitFile(filename: str):
    DIMENSIONS_SECTION = 6
    COORDINATES_SECTION = 9

    with open(filename, 'r') as f:
        lines = f.readlines()
    
    dims = lines[DIMENSIONS_SECTION].split()
    Nv = int(dims[0])
    Nk =  int(dims[1])
        
    vx = np.zeros(Nv, dtype=float)
    vy = np.zeros(Nv, dtype=float)
    for i in range(Nv):
        splitted_line = lines[i + COORDINATES_SECTION].split()
        vx[i] = float(splitted_line[1])
        vy[i] = float(splitted_line[2])
    
    EToV = np.zeros((Nk, N_FACES), dtype=int)
    for k in range(Nk):
        elements_section_begin = COORDINATES_SECTION + Nv + 2
        splitted_line = lines[k + elements_section_begin].split()
        EToV[k,0] = int(splitted_line[3]) - 1
        EToV[k,1] = int(splitted_line[4]) - 1
        EToV[k,2] = int(splitted_line[5]) - 1

    return Mesh2D(vx, vy, EToV)