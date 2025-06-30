import numpy as np
import matplotlib.pyplot as plt

from .dg3d_tools import *
from .mesh3d import Mesh3D
from ..integrators.LSERK4 import *
from ..spatialDiscretization import *


class Maxwell3D(SpatialDiscretization):
    def __init__(self, n_order: int, mesh: Mesh3D, fluxType="Upwind"):
        assert n_order > 0
        assert mesh.number_of_elements() > 0

        self.n_order = n_order
        self.n_fp = int((n_order + 1) * (n_order + 2) / 2)
        self.n_faces = 4
        
        self.mesh = mesh
        self.fluxType = fluxType

        self.epsilon = np.ones(mesh.number_of_elements())
        self.mu = np.ones(mesh.number_of_elements())

        r, s, t = xyz_to_rst(*set_nodes_in_equilateral_tetrahedron(n_order))
        self.V = vandermonde3D(n_order, r, s, t)
        self.mass_matrix = mass_matrix(n_order, r, s, t)
        self.Dr, self.Ds, self.Dt = differentiation_matrices_3d(n_order, r, s, t, self.V)
        self.x, self.y, self.z = nodes_coordinates(n_order, mesh)

        self.lift = lift3D(n_order)
        self.EToE, self.EToF = mesh.connectivityMatrices()

        self.rx, self.sx, self.tx, self.ry, self.sy, self.ty, self.rz, self.sz, self.tz, self.jacobian = geometricFactors3D(
            self.x, self.y, self.z, self.Dr, self.Ds, self.Dt)
        
        self.fmask, _, _, _, _ = buildFMask(n_order)

        self.nx, self.ny, self.nz, sJ = normals3D(self.x, self.y, self.z, self.Dr, self.Ds, self.Dt, n_order)
        self.f_scale = sJ/self.jacobian[self.fmask.ravel('F')]

        self.buildMaps()

    def buildMaps(self):
        '''
        function [vmapM, vmapP, vmapB, mapB] = BuildMaps3D
        % function [vmapM, vmapP, vmapB, mapB] = BuildMaps3D
        % Purpose: Connectivity and boundary tables for nodes given
        %           in the K#of elements, each with N+1 degrees of freedom.
        '''
        K = self.mesh.number_of_elements()
        Np = self.number_of_nodes_per_element()
        Nfaces = self.n_faces
        Nfp = self.n_fp

        # number volume nodes consecutively
        node_ids = np.reshape(np.arange(K*Np), [Np, K], 'F')
        vmapM = np.full([Nfp, Nfaces, K], 0)
        vmapP = np.full([Nfp, Nfaces, K], 0)

        # find index of face nodes with respect to volume node ordering
        for k1 in range(K):
            for f1 in range(Nfaces):
                vmapM[:, f1, k1] = node_ids[self.fmask[:, f1], k1]

        tmp = np.ones(Nfp)
        for k1 in range(K):
            for f1 in range(Nfaces):
                # find neighbor
                k2 = self.EToE[k1, f1]
                f2 = self.EToF[k1, f1]

                # find volume node numbers of left and right nodes
                vidM = vmapM[:, f1, k1]
                vidP = vmapM[:, f2, k2]
                
                xM = np.outer(self.x.ravel('F')[vidM], tmp)
                yM = np.outer(self.y.ravel('F')[vidM], tmp)
                zM = np.outer(self.z.ravel('F')[vidM], tmp)
                
                xP = np.outer(self.x.ravel('F')[vidP], tmp)
                yP = np.outer(self.y.ravel('F')[vidP], tmp)
                zP = np.outer(self.z.ravel('F')[vidP], tmp)

                # Compute distance matrix
                dist2 = (xM - xP.T)**2 + (yM - yP.T)**2 + (zM - zP.T)**2

                idM, idP = np.where(dist2 <= NODETOL)
                vmapP[idM, f1, k1] = vidP[idP]

        self.vmapM = vmapM.ravel('F')
        self.vmapP = vmapP.ravel('F')
        self.vmapB = self.vmapM[self.vmapP == self.vmapM]
        self.mapB = np.where(self.vmapP == self.vmapM)[0]

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

    def buildEvolutionOperator(self):
        Np = self.number_of_nodes_per_element()
        K = self.mesh.number_of_elements()
        Nfields = 6  # Ex, Ey, Ez, Hx, Hy, Hz
        N = Nfields * Np * K # total DOFs
        A = np.zeros((N,N))

        for i in range(N):
            fields = self.buildFields()
            node = i % Np
            elem = int(np.floor(i / Np)) % K
            field_block = i // (Np * K)

            # Set the appropriate field based on the block index
            if field_block == 0:  # Ex
                fields['Ex'][node, elem] = 1.0
            elif field_block == 1:  # Ey
                fields['Ey'][node, elem] = 1.0
            elif field_block == 2:  # Ez
                fields['Ez'][node, elem] = 1.0
            elif field_block == 3:  # Hx
                fields['Hx'][node, elem] = 1.0
            elif field_block == 4:  # Hy
                fields['Hy'][node, elem] = 1.0
            elif field_block == 5:  # Hz
                fields['Hz'][node, elem] = 1.0
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
                fieldsRHS['Hz'].reshape(Np*K, 1, order='F')
            ])

            # Fill the operator matrix with the computed values
            A[:,i] = q0[:,0]

        return A

    def buildFields(self):
        Ex = np.zeros([self.number_of_nodes_per_element(),
                       self.mesh.number_of_elements()])
        Ey = np.zeros(Ex.shape)
        Ez = np.zeros(Ex.shape)
        Hx = np.zeros(Ex.shape)
        Hy = np.zeros(Ex.shape)
        Hz = np.zeros(Ex.shape)

        return {'Hx': Hx, 'Hy': Hy, 'Hz': Hz,
                'Ex': Ex, 'Ey': Ey, 'Ez': Ez}

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

        dHx[self.mapB] = Hx.transpose().take(self.vmapB) - Hbcx
        dHy[self.mapB] = Hy.transpose().take(self.vmapB) - Hbcy
        dHz[self.mapB] = Hz.transpose().take(self.vmapB) - Hbcz
        
        dEx[self.mapB] = Ex.transpose().take(self.vmapB) - Ebcx
        dEy[self.mapB] = Ey.transpose().take(self.vmapB) - Ebcy
        dEz[self.mapB] = Ez.transpose().take(self.vmapB) - Ebcz

        dHx = dHx.reshape(self.n_fp*self.n_faces, self.mesh.number_of_elements(), order='F')
        dHy = dHy.reshape(self.n_fp*self.n_faces, self.mesh.number_of_elements(), order='F')
        dHz = dHz.reshape(self.n_fp*self.n_faces, self.mesh.number_of_elements(), order='F')
        
        dEx = dEx.reshape(self.n_fp*self.n_faces, self.mesh.number_of_elements(), order='F')
        dEy = dEy.reshape(self.n_fp*self.n_faces, self.mesh.number_of_elements(), order='F')
        dEz = dEz.reshape(self.n_fp*self.n_faces, self.mesh.number_of_elements(), order='F')

        return dHx, dHy, dHz, dEx, dEy, dEz

    def computeRHS(self, fields):
        Hx, Hy, Hz = fields['Hx'], fields['Hy'], fields['Hz']
        Ex, Ey, Ez = fields['Ex'], fields['Ey'], fields['Ez']

        flux_Hx, flux_Hy, flux_Hz, flux_Ex, flux_Ey, flux_Ez = self.computeFlux(Hx, Hy, Hz, Ex, Ey, Ez)

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

        return {'Hx': rhs_Hx, 'Hy': rhs_Hy, 'Hz': rhs_Hz,
                'Ex': rhs_Ex, 'Ey': rhs_Ey, 'Ez': rhs_Ez}
    