import numpy as np

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

    

