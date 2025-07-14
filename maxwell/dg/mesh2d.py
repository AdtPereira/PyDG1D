import numpy as np
import gmsh
import matplotlib.tri as mtri

N_FACES = 3

class Mesh2D:
    def __init__(self, vx, vy, EToV, boundary_label="PEC"):
        assert vx.shape == vy.shape
        assert np.max(np.max(EToV)) + 1  == vx.shape[0]
        
        self.dimension = 2
        self.vx = vx
        self.vy = vy
        self.EToV = EToV        
        self.boundary_label = boundary_label

    def number_of_vertices(self):
        return self.vx.shape[0]

    def number_of_elements(self):
        return self.EToV.shape[0]
    
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
            (self.EToV[:,[0,1]], 
             self.EToV[:,[1,2]], 
             self.EToV[:,[2,0]])
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
        self.physical_groups = None


    def __enter__(self):
        gmsh.initialize()
        gmsh.open(self.filename)
        return self


    def __exit__(self, exc_type, exc_val, exc_tb):
        gmsh.finalize()


    def read_global_mesh(self,  GroupDimTag: tuple[int, int] = (2, 201)) -> None:
        """
        Lê os dados da malha global do arquivo e os armazena nos atributos da instância.
        Este método preenche self.vx, self.vy e self.EToV para a malha inteira.
        """
        PhysicalGroups = self._getPhysicalNodesDict()
        nodesData = PhysicalGroups[GroupDimTag[1]]['nodes']
        self.vx, self.vy = self._extractVertices(nodesData)
        self.EToV = self._getEToV(GroupDimTag)
        self.physical_groups = gmsh.model.getPhysicalGroups()


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


    def _getPhysicalNodesDict(self):
        PhysicalGroups = {}
        for dim, GroupTag in gmsh.model.getPhysicalGroups():
            NodeTags, NodeCoords = gmsh.model.mesh.getNodesForPhysicalGroup(dim, GroupTag)
            NodeTags = [int(tag) for tag in NodeTags.astype(np.uint64)]

            PhysicalGroups[GroupTag] = {
                "name": gmsh.model.getPhysicalName(dim, GroupTag),
                "nodes": {
                    int(node): (
                        float(NodeCoords[3*i]),
                        float(NodeCoords[3*i + 1]),
                        float(NodeCoords[3*i + 2]),
                    ) for i, node in enumerate(NodeTags)
                }
            }

        return PhysicalGroups
    
    
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
        physical_group_data = {}
        # Coleta informações de todos os grupos de uma vez
        all_nodes_by_group = self._getPhysicalNodesDict()

        # Itera sobre a lista de grupos físicos, que já contém a dimensão e a tag
        for dim, tag in self.physical_groups:
            # Pula se, por algum motivo, o grupo não tiver nós associados
            if tag not in all_nodes_by_group:
                continue

            group_info = all_nodes_by_group[tag]
            group_name = group_info['name']
            nodes_data = group_info['nodes']

            # Extrai os vértices e a conectividade para este grupo específico
            # A função _getEToV funciona para qualquer dimensão
            EToV = self._getEToV((dim, tag))
            vx, vy = self._extractVertices(nodes_data)

            # Armazena os dados, incluindo a dimensão do grupo
            physical_group_data[group_name] = {
                'dim': dim,  # <-- PONTO-CHAVE: Armazena a dimensão do grupo
                'vx': vx,
                'vy': vy,
                'EToV': EToV # Nome mais genérico que 'EToV'
            }
        
        return physical_group_data
    

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