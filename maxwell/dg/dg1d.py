import numpy as np
import math
from .dg1d_tools import *
from ..spatialDiscretization import *
from .mesh1d import Mesh1D


class DG1D(SpatialDiscretization):
    def __init__(self, n_order: int, mesh: Mesh1D, fluxType="Upwind",epsilon=None,sigma=None):
        SpatialDiscretization.__init__(self, mesh)
        
        assert n_order > 0
        self.n_order = n_order        
        assert fluxType == "Upwind" or fluxType == "Centered"
        self.fluxType = fluxType       
        self.n_faces = 2
        self.n_fp = 1   

        alpha = 0
        beta = 0

        # Epsilon implementation in 1D
        if epsilon is None:
            self.epsilon = np.ones(mesh.number_of_elements())
        elif len(epsilon) != mesh.number_of_elements():
            raise ValueError("The dimensions of the permittivity vector must align with the number of elements in the mesh.")
        else:          
            self.epsilon = np.array(epsilon)


        # Sigma implementation necessary for J for 1D
        if sigma is None:
            self.sigma = np.zeros(mesh.number_of_elements())
        elif len(sigma) != mesh.number_of_elements():
            raise ValueError("The dimensions of the charge density vector must align with the number of elements in the mesh.")
        else:          
            self.sigma = np.array(sigma)       
        
        
        self.mu = np.ones(mesh.number_of_elements())
        self.x = nodes_coordinates(n_order, mesh.EToV, mesh.vx)
        self.nx = normals(mesh.number_of_elements())

        etoe, etof = connect(mesh.EToV)
        self.vmap_m, self.vmap_p, self.vmap_b, self.map_b = build_maps(
            n_order, self.x, etoe, etof)

        self.r = jacobiGL(alpha, beta, n_order)
        self.fmask, self.fmask_1, self.fmask_2 = buildFMask(self.r)

        self.mass = mass_matrix(n_order, self.r)
        self.lift = surface_integral_dg(n_order, self.r)
        self.diff_matrix = differentiation_matrix(n_order, self.r)

        self.rx, self.jacobian = geometric_factors(self.x, self.diff_matrix)
        self.f_scale = 1/self.jacobian[self.fmask]

        K = self.mesh.number_of_elements()
        Z_imp = self.get_impedance()

        self.Z_imp_m = Z_imp.transpose().take(self.vmap_m)
        self.Z_imp_p = Z_imp.transpose().take(self.vmap_p)
        self.Z_imp_m = self.Z_imp_m.reshape(
            self.n_fp*self.n_faces, K, order='F')
        self.Z_imp_p = self.Z_imp_p.reshape(
            self.n_fp*self.n_faces, K, order='F')

        self.Y_imp_m = 1.0 / self.Z_imp_m
        self.Y_imp_p = 1.0 / self.Z_imp_p

        self.Z_imp_sum = self.Z_imp_m + self.Z_imp_p
        self.Y_imp_sum = self.Y_imp_m + self.Y_imp_p

    def number_of_nodes_per_element(self):
        return self.n_order + 1

    def get_nodes(self):
        return set_nodes(self.n_order, self.mesh.vx[self.mesh.EToV])

    def get_minimum_node_distance(self):
        return min(np.abs(self.x[0, :] - self.x[1, :]))

    def buildFields(self):
        E = np.zeros([self.number_of_nodes_per_element(),
                      self.mesh.number_of_elements()])
        H = np.zeros(E.shape)


        return {"E": E, "H": H}

    def get_impedance(self):
        Z_imp = np.zeros(self.x.shape)
        for i in range(Z_imp.shape[1]):
            Z_imp[:, i] = np.sqrt(self.mu[i] / self.epsilon[i])

        return Z_imp

    def fieldsOnBoundaryConditions(self, E, H):
        bcType = self.mesh.boundary_label

        for bdr, label in self.mesh.boundary_label.items():
            if bdr == "LEFT" or bdr == "RIGHT":
                if label == "PEC":
                    Ebc = - E.transpose().take(self.vmap_b)
                    Hbc = H.transpose().take(self.vmap_b)
                elif label == "PMC":
                    Hbc = - H.transpose().take(self.vmap_b)
                    Ebc = E.transpose().take(self.vmap_b)
                elif label == "SMA":
                    Hbc = H.transpose().take(self.vmap_b) * 0.0
                    Ebc = E.transpose().take(self.vmap_b) * 0.0
                elif label == "Periodic":
                    Ebc = E.transpose().take(self.vmap_b[::-1])
                    Hbc = H.transpose().take(self.vmap_b[::-1])
                else:
                    raise ValueError("Invalid boundary label.")
                return Ebc, Hbc

    def computeFluxE(self, E, H):
        dE, dH = self.computeJumps(E, H)

        if self.fluxType == "Upwind":
            flux_E = 1/self.Z_imp_sum*(self.nx*self.Z_imp_p*dH-dE)
        elif self.fluxType == "Centered":
            flux_E = 1/self.Z_imp_sum*(self.nx*self.Z_imp_p*dH)
        else:
            raise ValueError("Invalid fluxType label")
        return flux_E

    def computeFluxH(self, E, H):
        dE, dH = self.computeJumps(E, H)

        if self.fluxType == "Upwind":
            flux_H = 1/self.Y_imp_sum*(self.nx*self.Y_imp_p*dE-dH)
        elif self.fluxType == "Centered":
            flux_H = 1/self.Y_imp_sum*(self.nx*self.Y_imp_p*dE)
        else:
            raise ValueError("Invalid fluxType label")
        return flux_H

    def computeFlux(self, E, H):
        dE, dH = self.computeJumps(E, H)

        if self.fluxType == "Upwind":
            flux_E = 1/self.Z_imp_sum*(self.nx*self.Z_imp_p*dH-dE)
            flux_H = 1/self.Y_imp_sum*(self.nx*self.Y_imp_p*dE-dH)
        elif self.fluxType == "Centered":
            flux_E = 1/self.Z_imp_sum*(self.nx*self.Z_imp_p*dH)
            flux_H = 1/self.Y_imp_sum*(self.nx*self.Y_imp_p*dE)
        else:
            raise ValueError("Invalid fluxType label")
        return flux_E, flux_H

    def computeJumps(self, E, H):
        Ebc, Hbc = self.fieldsOnBoundaryConditions(E, H)
        dE = E.transpose().take(self.vmap_m) - E.transpose().take(self.vmap_p)
        dH = H.transpose().take(self.vmap_m) - H.transpose().take(self.vmap_p)
        dE[self.map_b] = E.transpose().take(self.vmap_b)-Ebc
        dH[self.map_b] = H.transpose().take(self.vmap_b)-Hbc
        dE = dE.reshape(self.n_fp*self.n_faces,
                        self.mesh.number_of_elements(), order='F')
        dH = dH.reshape(self.n_fp*self.n_faces,
                        self.mesh.number_of_elements(), order='F')
        return dE, dH

    def computeRHSE(self, fields):
        E = fields['E']
        H = fields['H']
        J = np.zeros((self.number_of_nodes_per_element(), self.mesh.number_of_elements()))

        J[:, :] = E * self.sigma

        flux_E = self.computeFluxE(E, H)
        rhs_drH = np.matmul(self.diff_matrix, H)
        rhsE = 1/self.epsilon * \
            (np.multiply(-1*self.rx, rhs_drH) +
             np.matmul(self.lift, self.f_scale * flux_E)-J) 

        return rhsE

    def computeRHSH(self, fields):
        E = fields['E']
        H = fields['H']

        flux_H = self.computeFluxH(E, H)
        rhs_drE = np.matmul(self.diff_matrix, E)
        rhsH = 1/self.mu * (np.multiply(-1*self.rx, rhs_drE) +
                            np.matmul(self.lift, self.f_scale * flux_H))
        
        return rhsH

    def computeRHS(self, fields):
        rhsE = self.computeRHSE(fields)
        rhsH = self.computeRHSH(fields)

        return {'E': rhsE, 'H': rhsH}

    def convertToVector(self, fields):
        return np.concatenate((
            fields['E'].ravel(order='F'),
            fields['H'].ravel(order='F')
        ))

    def copyVectorToFields(self, vec, fields):
        Np = self.number_of_nodes_per_element()
        K = self.mesh.number_of_elements()
        fields['E'][:, :] = vec[:(vec.size//2)].reshape(Np, K, order='F')
        fields['H'][:, :] = vec[(vec.size//2):].reshape(Np, K, order='F')

    def setFieldWithIndex(self, fields, i, val):
        Np = self.number_of_nodes_per_element()
        node = i % Np
        elem = int(np.floor(i / Np)) % self.mesh.number_of_elements()
        if i < self.number_of_unknowns()/2:
            fields['E'][node, elem] = val
        else:
            fields['H'][node, elem] = val
        return fields

    def buildEvolutionOperator(self):
        Np = self.number_of_nodes_per_element()
        K = self.mesh.number_of_elements()
        N = self.number_of_unknowns()
        A = np.zeros((N, N))
        for i in range(N):
            fields = self.buildFields()
            self.setFieldWithIndex(fields, i, 1.0)
            fieldsRHS = self.computeRHS(fields)
            q0 = np.vstack([
                fieldsRHS['E'].reshape(Np*K, 1, order='F'),
                fieldsRHS['H'].reshape(Np*K, 1, order='F')
            ])
            A[:, i] = q0[:, 0]
        return A
    
    def reorder_array(self, A, ordering):
        # Assumes that the original array contains all DoF ordered as:
        # [ E_0, ..., E_{K-1}, H_0, ..., H_{K-1} ]
        N = A.shape[0]
        K = self.mesh.number_of_elements()
        Np = self.number_of_nodes_per_element()
        new_order = np.arange(N, dtype=int)
        if ordering == 'byElements':
            for i in range(N):
                node = i % Np
                elem = int(np.floor(i / Np)) % K
                if i < N/2:
                    new_order[2*elem*Np+node] = i
                else:
                    new_order[2*elem*Np+Np+node] = i
        if ordering == 'interleaved':
            for i in range(N):
                node = i % Np
                elem = int(np.floor(i / Np)) % K
                if i < N/2:
                    new_order[2*elem*Np+node*2] = i
                else:
                    new_order[2*elem*Np+node*2+1] = i
        if (len(A.shape) == 1):
            A1 = [A[i] for i in new_order]
        elif (len(A.shape) == 2):
            A1 = [[A[i][j] for j in new_order] for i in new_order]
        return np.array(A1)

    def buildGlobalMassMatrix(self):
        Np = self.number_of_nodes_per_element()
        K = self.mesh.number_of_elements()
        N = 2 * Np * K
        M = np.zeros((N, N))
        for k in range(K):
            ini = k*Np
            end = (k+1)*Np
            M[ini:end, ini:end] = self.mass * self.jacobian[0, k]
        M[int(N/2):, int(N/2):] = M[:int(N/2), :int(N/2)]

        return M

    def getEnergy(self, field):
        '''
        Gets energy stored in field by computing
            field^T * MassMatrix * field * Jacobian.
        for each element and then the sum.
        '''
        Np = self.number_of_nodes_per_element()
        K = self.mesh.number_of_elements()
        assert field.shape == (Np, K)
        energy = 0.0
        for k in range(K):
            energy += np.inner(
                field[:, k].dot(self.mass),
                field[:, k]*self.jacobian[:, k]
            )

        return energy

    def getEnergyPerElement(self, field):
        """
        Retorna a energia armazenada em cada elemento separadamente:
            E_k = field^T * MassMatrix * field * Jacobian_k
        """
        Np = self.number_of_nodes_per_element()
        K = self.mesh.number_of_elements()
        assert field.shape == (Np, K)
        energy_per_element = np.zeros(K)

        for k in range(K):
            v = field[:, k]
            energy_per_element[k] = v @ (self.mass @ v) * self.jacobian[0, k]

        return energy_per_element