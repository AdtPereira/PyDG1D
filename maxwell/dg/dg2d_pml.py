import numpy as np
import matplotlib.pyplot as plt

from .dg2d_tools import *
from .dg1d_tools import jacobi_gauss
from .mesh2d import Mesh2D
from ..integrators.LSERK4 import *
from ..spatialDiscretization import *


class Maxwell2D(SpatialDiscretization):
    def __init__(self, n_order: int, mesh: Mesh2D, fluxType="Upwind", epsilon=None, sigma=None, pml_design=None):
        assert n_order > 0
        assert mesh.number_of_elements() > 0

        self.n_order = n_order
        self.n_fp = n_order + 1
        self.n_faces = 3
        
        self.mesh = mesh
        self.fluxType = fluxType

        # Epsilon implementation
        if epsilon is None:
            self.epsilon = np.ones(mesh.number_of_elements())
        elif len(epsilon) != mesh.number_of_elements():
            raise ValueError("The dimensions of the permittivity vector must align with the number of elements in the mesh.")
        else:          
            self.epsilon = np.array(epsilon)
        
        # Mu implementation
        self.mu = np.ones(mesh.number_of_elements())

        # Sigma implementation necessary for J
        if sigma is None:
            self.sigma = np.zeros(mesh.number_of_elements())
        elif len(sigma) != mesh.number_of_elements():
            raise ValueError("The dimensions of the charge density vector must align with the number of elements in the mesh.")
        else:          
            self.sigma = np.array(sigma) 

        r, s = xy_to_rs(*set_nodes_in_equilateral_triangle(n_order))
        self.mass = mass_matrix(n_order, r, s)
        self.Dr, self.Ds = derivateMatrix(n_order, r, s)
        self.x, self.y = nodes_coordinates(n_order, mesh)

        self.lift = lift(n_order)

        eToE, eToF = mesh.connectivityMatrices()
        va = self.mesh.EToV[:, 0]
        vb = self.mesh.EToV[:, 1]
        vc = self.mesh.EToV[:, 2]

        self.rx, self.sx, self.ry, self.sy, self.jacobian = geometricFactors(
            self.x, self.y, self.Dr, self.Ds)
        
        self.nx, self.ny, sJ = normals(self.x, self.y, self.Dr, self.Ds, n_order)
        
        fmask, _, _, _ = buildFMask(n_order)
        self.f_scale = sJ/self.jacobian[fmask.ravel('F')]

        self.buildMaps()

        # Build PML fields
        if pml_design is not None:
            self.sigma_x, self.sigma_y, self.dsigma_dx, self.dsigma_dy = self.evaluate_sigma_fields(
                self.x, self.y, pml_design)
        else:
            self.sigma_x = np.zeros_like(self.x)
            self.sigma_y = np.zeros_like(self.y)
            self.dsigma_dx = np.zeros_like(self.x)
            self.dsigma_dy = np.zeros_like(self.y)
            

    def evaluate_sigma_fields(self, x, y, pml_design):
        """
        Avalia os campos σx(x) e σy(y) para as regiões de PML.
        Parâmetros:
        - x, y: arrays de coordenadas onde avaliar os campos.
        - x_limit, y_limit: limites do domínio físico onde σ = 0.
        - sigma0_x, sigma0_y: coeficientes máximos de amortecimento.
        - p: expoente do perfil PML.
        Retorna:
        - sigma_x, sigma_y: arrays com os valores ponto a ponto para cada elemento.
        - dsigma_dx, dsigma_dy: derivadas dos campos σx e σy.
        """
        # σx(x), σy(y) e suas derivadas
        sigma_x = np.zeros_like(x)
        sigma_y = np.zeros_like(y)
        dsigma_dx = np.zeros_like(x)
        dsigma_dy = np.zeros_like(y)

        # Valor inicial e perfil PML
        sigma0_x = -np.log(pml_design['R'])
        sigma0_y = sigma0_x
        p = pml_design['pml_order']

        # Verifica se os limites são válidos
        x0, y0 = pml_design['x0'], pml_design['x0']
        if x0 <= 0 or y0 <= 0:
            raise ValueError("Os limites x_limit e y_limit devem ser positivos.")

        # σ^x(x)
        x_r = x >= x0
        x_l = x <= -x0
        sigma_x[x_r] = sigma0_x * (x[x_r] - x0)**p
        sigma_x[x_l]  = sigma0_x * (x[x_l] + x0)**p
        dsigma_dx[x_r] = p * sigma0_x * (x[x_r] - x0)**(p - 1)
        dsigma_dx[x_l]  = p * sigma0_x * (x[x_l] + x0)**(p - 1)

        # σ^y(y)
        y_t = y >= y0
        y_b = y <= -y0
        sigma_y[y_t]    = sigma0_y * (y[y_t] - y0)**p
        sigma_y[y_b] = sigma0_y * (y[y_b] + y0)**p
        dsigma_dy[y_t]    = p * sigma0_y * (y[y_t] - y0)**(p - 1)
        dsigma_dy[y_b] = p * sigma0_y * (y[y_b] + y0)**(p - 1)

        return sigma_x, sigma_y, dsigma_dx, dsigma_dy
    
    def buildMaps(self):
        '''
        function [mapM, mapP, vmapM, vmapP, vmapB, mapB] = BuildMaps2D
        Purpose: Connectivity and boundary tables in the K # of Np elements        
        '''
        N = self.n_order
        msh = self.mesh
        k_elem = self.mesh.number_of_elements()
        n_p = self.number_of_nodes_per_element()
        n_faces = 3
        n_fp = N+1

        # mask defined in globals
        Fmask, _, _, _ = buildFMask(N)

        # number volume nodes consecutively
        node_ids = np.reshape(np.arange(k_elem*n_p), [n_p, k_elem], 'F')
        vmapM = np.full([n_fp, n_faces, k_elem], 0)
        vmapP = np.full([n_fp, n_faces, k_elem], 0)
        mapM = np.arange(k_elem*n_fp*n_faces)
        mapP = np.reshape(mapM, (n_fp, n_faces, k_elem))

        # find index of face nodes with respect to volume node ordering
        for k1 in range(k_elem):
            for f1 in range(n_faces):
                vmapM[:, f1, k1] = node_ids[Fmask[:, f1], k1]

        one = np.ones(n_fp)
        EToE, EToF = msh.connectivityMatrices()
        for k1 in range(k_elem):
            for f1 in range(n_faces):
                # find neighbor
                k2 = EToE[k1, f1]
                f2 = EToF[k1, f1]

                # reference length of edge
                v1 = msh.EToV[k1, f1]
                v2 = msh.EToV[k1, np.mod(f1+1, n_faces)]
                refd = np.sqrt(
                    (msh.vx[v1]-msh.vx[v2])**2 + (msh.vy[v1]-msh.vy[v2])**2
                )

                # find find volume node numbers of left and right nodes
                vidM = vmapM[:, f1, k1]
                vidP = vmapM[:, f2, k2]
                x1 = np.outer(self.x.ravel('F')[vidM], one)
                y1 = np.outer(self.y.ravel('F')[vidM], one)
                x2 = np.outer(self.x.ravel('F')[vidP], one)
                y2 = np.outer(self.y.ravel('F')[vidP], one)

                # Compute distance matrix
                distance = np.sqrt(
                    np.abs((x1 - x2.transpose())**2 + (y1-y2.transpose())**2))
                idM, idP = np.where(distance <= NODETOL*refd)
                vmapP[idM, f1, k1] = vidP[idP]
                mapP[idM, f1, k1] = idP + (f2-1)*n_fp+(k2-1)*n_faces*n_fp

        vmapM = vmapM.ravel('F')
        vmapP = vmapP.ravel('F')
        vmapB = vmapM[vmapP == vmapM]
        mapB = np.where(vmapP == vmapM)[0]

        self.vmapM = vmapM
        self.vmapP = vmapP
        self.vmapB = vmapB
        self.mapB = mapB

    def get_minimum_node_distance(self):
        points, _ = jacobi_gauss(0, 0, self.n_order)
        return abs(points[0]-points[1])
    
    def get_dt_scale(self):

        r, s = xy_to_rs(*set_nodes_in_equilateral_triangle(self.n_order))
        vmask1 = np.where(np.abs(s+r+2) < NODETOL)[0]
        vmask2 = np.where(np.abs(r-1) < NODETOL)[0]
        vmask3 = np.where(np.abs(s-1) < NODETOL)[0]
        vmask  = np.array([vmask1, vmask2, vmask3]).transpose()

        vx = self.x[np.squeeze(vmask.reshape(-1, 1)), :]
        vy = self.y[np.squeeze(vmask.reshape(-1, 1)), :]

        len1 = np.sqrt((vx[0,:]-vx[1,:])**2+(vy[0,:]-vy[1,:])**2)
        len2 = np.sqrt((vx[1,:]-vx[2,:])**2+(vy[1,:]-vy[2,:])**2)
        len3 = np.sqrt((vx[2,:]-vx[0,:])**2+(vy[2,:]-vy[0,:])**2)
        sper = (len1 + len2 + len3)/2.0
        area = np.sqrt(sper*(sper-len1)*(sper-len2)*(sper-len3))

        dtscale = area/sper

        return dtscale

    def get_mesh(self):
        return self.mesh

    def number_of_nodes_per_element(self):
        return int((self.n_order + 1) * (self.n_order + 2) / 2)
    
    def buildEvolutionOperator(self):
        Np = self.number_of_nodes_per_element()
        K = self.mesh.number_of_elements()
        N_fields = 7
        N = N_fields * Np * K        
        A = np.zeros((N,N))

        for i in range(N):
            fields = self.buildFields()
            node = i % Np
            elem = int(np.floor(i / Np)) % K

            block = i // (Np * K) 
            field_names = ['Ez', 'Hx', 'Hy', 'Px', 'Py', 'Qx', 'Qy']
            fields[field_names[block]][node, elem] = 1.0

            # RHS completo com termos PML
            rhs = self.computeRHS(fields)

            # Empilhar todos os campos na mesma ordem
            q_rhs = np.vstack([
                rhs['Ez'].reshape(Np * K, 1, order='F'),
                rhs['Hx'].reshape(Np * K, 1, order='F'),
                rhs['Hy'].reshape(Np * K, 1, order='F'),
                rhs['Px'].reshape(Np * K, 1, order='F'),
                rhs['Py'].reshape(Np * K, 1, order='F'),
                rhs['Qx'].reshape(Np * K, 1, order='F'),
                rhs['Qy'].reshape(Np * K, 1, order='F')
            ])

            A[:, i] = q_rhs[:, 0]

        return A

    def buildStiffnessEvolutionOperator(self):
        Np = self.number_of_nodes_per_element()
        K = self.mesh.number_of_elements()
        N = 3 * Np * K
        A = np.zeros((N,N))
        for i in range(N):
            fields = self.buildFields()
            node = i % Np
            elem = int(np.floor(i / Np)) % K
            if i < N/3:
                fields['Ez'][node,elem] = 1.0
            elif i >= N/3 and i < 2*N/3:
                fields['Hx'][node,elem] = 1.0
            else:
                fields['Hy'][node,elem] = 1.0
            fieldsStiffness = self.computeRHSStiffness(fields)
            q0 = np.vstack([
                fieldsStiffness['Ez'].reshape(Np*K,1,order='F'),
                fieldsStiffness['Hx'].reshape(Np*K,1,order='F'),
                fieldsStiffness['Hy'].reshape(Np*K,1,order='F')
            ])
            A[:,i] = q0[:,0]

        return A

    def buildZeroNormalEvolutionOperator(self):
        Np = self.number_of_nodes_per_element()
        K = self.mesh.number_of_elements()
        N = 3 * Np * K
        A = np.zeros((N,N))
        for i in range(N):
            fields = self.buildFields()
            node = i % Np
            elem = int(np.floor(i / Np)) % K
            if i < N/3:
                fields['Ez'][node,elem] = 1.0
            elif i >= N/3 and i < 2*N/3:
                fields['Hx'][node,elem] = 1.0
            else:
                fields['Hy'][node,elem] = 1.0
            fieldsZeroNormal = self.computeRHSZeroNormal(fields)
            q0 = np.vstack([
                fieldsZeroNormal['Ez'].reshape(Np*K,1,order='F'),
                fieldsZeroNormal['Hx'].reshape(Np*K,1,order='F'),
                fieldsZeroNormal['Hy'].reshape(Np*K,1,order='F')
            ])
            A[:,i] = q0[:,0]

        return A
    
    def buildOneNormalEvolutionOperator(self):
        Np = self.number_of_nodes_per_element()
        K = self.mesh.number_of_elements()
        N = 3 * Np * K
        A = np.zeros((N,N))
        for i in range(N):
            fields = self.buildFields()
            node = i % Np
            elem = int(np.floor(i / Np)) % K
            if i < N/3:
                fields['Ez'][node,elem] = 1.0
            elif i >= N/3 and i < 2*N/3:
                fields['Hx'][node,elem] = 1.0
            else:
                fields['Hy'][node,elem] = 1.0
            fieldsOneNormal = self.computeRHSOneNormal(fields)
            q0 = np.vstack([
                fieldsOneNormal['Ez'].reshape(Np*K,1,order='F'),
                fieldsOneNormal['Hx'].reshape(Np*K,1,order='F'),
                fieldsOneNormal['Hy'].reshape(Np*K,1,order='F')
            ])
            A[:,i] = q0[:,0]

        return A
    
    def buildTwoNormalEvolutionOperator(self):
        Np = self.number_of_nodes_per_element()
        K = self.mesh.number_of_elements()
        N = 3 * Np * K
        A = np.zeros((N,N))
        for i in range(N):
            fields = self.buildFields()
            node = i % Np
            elem = int(np.floor(i / Np)) % K
            if i < N/3:
                fields['Ez'][node,elem] = 1.0
            elif i >= N/3 and i < 2*N/3:
                fields['Hx'][node,elem] = 1.0
            else:
                fields['Hy'][node,elem] = 1.0
            fieldsTwoNormal = self.computeRHSTwoNormal(fields)
            q0 = np.vstack([
                fieldsTwoNormal['Ez'].reshape(Np*K,1,order='F'),
                fieldsTwoNormal['Hx'].reshape(Np*K,1,order='F'),
                fieldsTwoNormal['Hy'].reshape(Np*K,1,order='F')
            ])
            A[:,i] = q0[:,0]

        return A

    def buildFields(self):
        Hx = np.zeros([self.number_of_nodes_per_element(),
                       self.mesh.number_of_elements()])
        Hy = np.zeros(Hx.shape)
        Ez = np.zeros(Hx.shape)
        Px = np.zeros(Hx.shape)
        Py = np.zeros(Hx.shape)
        Qx = np.zeros(Hx.shape)
        Qy = np.zeros(Hx.shape)

        return {'Hx': Hx, 'Hy': Hy, 'Ez': Ez, 'Px': Px, 'Py': Py, 'Qx': Qx, 'Qy': Qy}

    def computeZeroNormalFlux(self, dEz):

        flux_Hx_Zero_Normal = 0
        flux_Hy_Zero_Normal = 0
        flux_Ez_Zero_Normal = 0

        if self.fluxType == "Upwind":
            flux_Ez_Zero_Normal -= dEz
        elif self.fluxType == "Centered":
            pass
        else:
            raise ValueError("Invalid flux type.")
        
        return flux_Hx_Zero_Normal, flux_Hy_Zero_Normal, flux_Ez_Zero_Normal
    
    def computeOneNormalFlux(self, dHx, dHy, dEz):

        flux_Hx_One_Normal =  self.ny * dEz
        flux_Hy_One_Normal = -self.nx * dEz
        flux_Ez_One_Normal = -self.nx * dHy + self.ny * dHx

        return flux_Hx_One_Normal, flux_Hy_One_Normal, flux_Ez_One_Normal
    
    def computeTwoNormalFlux(self, dHx, dHy):

        flux_Hx_Two_Normal = 0
        flux_Hy_Two_Normal = 0
        flux_Ez_Two_Normal = 0

        if self.fluxType == "Upwind":
            ndotdH = self.nx * dHx + self.ny * dHy
            flux_Hx_Two_Normal += ndotdH * self.nx - dHx
            flux_Hy_Two_Normal += ndotdH * self.ny - dHy
        elif self.fluxType == "Centered":
            pass
        else:
            raise ValueError("Invalid flux type.")
        
        return flux_Hx_Two_Normal, flux_Hy_Two_Normal, flux_Ez_Two_Normal

    def computeFlux(self, Hx, Hy, Ez):

        dHx, dHy, dEz = self.computeJumps(Hx, Hy, Ez)
        f_Hx_zero, f_Hy_zero, f_Ez_zero = self.computeZeroNormalFlux(dEz)
        f_Hx_one,   f_Hy_one,  f_Ez_one = self.computeOneNormalFlux(dHx, dHy, dEz)
        f_Hx_two,   f_Hy_two,  f_Ez_two = self.computeTwoNormalFlux(dHx, dHy)
        flux_Hx =  f_Hx_zero + f_Hx_one + f_Hx_two
        flux_Hy =  f_Hy_zero + f_Hy_one + f_Hy_two
        flux_Ez =  f_Ez_zero + f_Ez_one + f_Ez_two

        return flux_Hx, flux_Hy, flux_Ez

    def fieldsOnBoundaryConditions(self, Hx, Hy, Ez):

        bcType = self.mesh.boundary_label
        if bcType == "PEC":
            Hbcx = Hx.transpose().take(self.vmapB)
            Hbcy = Hy.transpose().take(self.vmapB)
            Ebcz = - Ez.transpose().take(self.vmapB)
        elif bcType == "PMC":
            Hbcx = - Hx.transpose().take(self.vmapB)
            Hbcy = - Hy.transpose().take(self.vmapB)
            Ebcz = Ez.transpose().take(self.vmapB)
        elif bcType == "SMA":
            Hbcx = Hx.transpose().take(self.vmapB) * 0.0
            Hbcy = Hx.transpose().take(self.vmapB) * 0.0
            Ebcz = Ez.transpose().take(self.vmapB) * 0.0
        elif bcType == "Periodic":
            Hbcx = Hx.transpose().take(self.vmapB[::-1])
            Hbcy = Hy.transpose().take(self.vmapB[::-1])
            Ebcz = Ez.transpose().take(self.vmapB[::-1])
        else:
            raise ValueError("Invalid boundary label.")
        return Hbcx, Hbcy, Ebcz

    def computeJumps(self, Hx, Hy, Ez):
        Hbcx, Hbcy, Ebcz = self.fieldsOnBoundaryConditions(Hx, Hy, Ez)
        dHx = Hx.transpose().take(self.vmapM) - Hx.transpose().take(self.vmapP)
        dHy = Hy.transpose().take(self.vmapM) - Hy.transpose().take(self.vmapP)
        dEz = Ez.transpose().take(self.vmapM) - Ez.transpose().take(self.vmapP)

        dHx[self.mapB] = Hx.transpose().take(self.vmapB) - Hbcx
        dHy[self.mapB] = Hy.transpose().take(self.vmapB) - Hbcy
        dEz[self.mapB] = Ez.transpose().take(self.vmapB) - Ebcz

        dHx = dHx.reshape(self.n_fp*self.n_faces,
                          self.mesh.number_of_elements(), order='F')
        dHy = dHy.reshape(self.n_fp*self.n_faces,
                          self.mesh.number_of_elements(), order='F')
        dEz = dEz.reshape(self.n_fp*self.n_faces,
                          self.mesh.number_of_elements(), order='F')

        return dHx, dHy, dEz

    def computeRHS(self, fields):
        Hx, Hy, Ez = fields['Hx'], fields['Hy'], fields['Ez']
        Px, Py = fields['Px'], fields['Py']
        Qx, Qy = fields['Qx'], fields['Qy']
        
        # Compute the RHS for the time evolution
        flux_Hx, flux_Hy, flux_Ez = self.computeFlux(Hx, Hy, Ez)

        rhs_Ezx, rhs_Ezy = grad(self.Dr, self.Ds, Ez, self.rx, self.sx, self.ry, self.sy)
        rhs_CuHz = curl(self.Dr, self.Ds, Hx, Hy, self.rx, self.sx, self.ry, self.sy)

        # -------- Current Density J ----------
        Jz = np.zeros((self.number_of_nodes_per_element(), self.mesh.number_of_elements()))
        Jz[:, :] = Ez * self.sigma

        # -------- Termos extras da formulação PML  ----------
        # -------- missing material epsilon/mu      ----------
        rhs_Hx = -rhs_Ezy  - self.sigma_y * (2 * Hx + Py) + np.matmul(self.lift, self.f_scale * flux_Hx)/2.0
        rhs_Hy =  rhs_Ezx  - self.sigma_x * (2 * Hy + Px) + np.matmul(self.lift, self.f_scale * flux_Hy)/2.0
        rhs_Ez =  rhs_CuHz - self.dsigma_dx * Qx + self.dsigma_dy * Qy + np.matmul(self.lift, self.f_scale * flux_Ez)/2.0 - Jz
        
        rhs_Px =  self.sigma_x * Hy
        rhs_Py =  self.sigma_y * Hx
        rhs_Qx = -self.sigma_x * Qx - Hy
        rhs_Qy = -self.sigma_y * Qy - Hx

        return {'Hx': rhs_Hx, 'Hy': rhs_Hy, 'Ez': rhs_Ez, 'Px': rhs_Px, 'Py': rhs_Py, 'Qx': rhs_Qx, 'Qy': rhs_Qy}

    def computeRHSStiffness(self, fields):
        # Compute the RHS for the stiffness matrix
        Hx, Hy, Ez = fields['Hx'], fields['Hy'], fields['Ez']       

        rhs_Ezx, rhs_Ezy = grad(self.Dr, self.Ds, Ez, self.rx, self.sx, self.ry, self.sy)
        rhs_CuHz = curl(self.Dr, self.Ds, Hx, Hy, self.rx, self.sx, self.ry, self.sy)

        # -------- Termos extras da formulação PML ----------
        rhs_Hx_stiffness = -rhs_Ezy
        rhs_Hy_stiffness =  rhs_Ezx
        rhs_Ez_stiffness =  rhs_CuHz 

        return {'Hx': rhs_Hx_stiffness, 'Hy': rhs_Hy_stiffness, 'Ez': rhs_Ez_stiffness}
    
    def computeRHSZeroNormal(self, fields):
        Hx = fields['Hx']
        Hy = fields['Hy']
        Ez = fields['Ez']

        _ , _ , dEz = self.computeJumps(Hx, Hy, Ez)
        flux_Hx_zero_normal, flux_Hy_zero_normal, flux_Ez_zero_normal = self.computeZeroNormalFlux(dEz)

        rhs_Hx_zero_normal = np.matmul(self.lift, self.f_scale * flux_Hx_zero_normal)/2.0
        rhs_Hy_zero_normal = np.matmul(self.lift, self.f_scale * flux_Hy_zero_normal)/2.0
        rhs_Ez_zero_normal = np.matmul(self.lift, self.f_scale * flux_Ez_zero_normal)/2.0

        return {'Hx': rhs_Hx_zero_normal, 'Hy': rhs_Hy_zero_normal, 'Ez': rhs_Ez_zero_normal}
    
    def computeRHSOneNormal(self, fields):
        Hx = fields['Hx']
        Hy = fields['Hy']
        Ez = fields['Ez']

        dHx , dHy , dEz = self.computeJumps(Hx, Hy, Ez)
        flux_Hx_one_normal, flux_Hy_one_normal, flux_Ez_one_normal = self.computeOneNormalFlux(dHx , dHy , dEz)

        rhs_Hx_one_normal = np.matmul(self.lift, self.f_scale * flux_Hx_one_normal)/2.0
        rhs_Hy_one_normal = np.matmul(self.lift, self.f_scale * flux_Hy_one_normal)/2.0
        rhs_Ez_one_normal = np.matmul(self.lift, self.f_scale * flux_Ez_one_normal)/2.0

        return {'Hx': rhs_Hx_one_normal, 'Hy': rhs_Hy_one_normal, 'Ez': rhs_Ez_one_normal}
    
    def computeRHSTwoNormal(self, fields):
        Hx = fields['Hx']
        Hy = fields['Hy']
        Ez = fields['Ez']

        dHx , dHy , _ = self.computeJumps(Hx, Hy, Ez)
        flux_Hx_two_normal, flux_Hy_two_normal, flux_Ez_two_normal = self.computeTwoNormalFlux(dHx, dHy)

        rhs_Hx_two_normal = np.matmul(self.lift, self.f_scale * flux_Hx_two_normal)/2.0
        rhs_Hy_two_normal = np.matmul(self.lift, self.f_scale * flux_Hy_two_normal)/2.0
        rhs_Ez_two_normal = np.matmul(self.lift, self.f_scale * flux_Ez_two_normal)/2.0

        return {'Hx': rhs_Hx_two_normal, 'Hy': rhs_Hy_two_normal, 'Ez': rhs_Ez_two_normal}

    def plot_field(self, Nout, field):
        # Build equally spaced grid on reference triangle
        Npout = int((Nout+1)*(Nout+2)/2)
        rout = np.zeros((Npout))
        sout = np.zeros((Npout))
        counter = np.zeros((Nout+1, Nout+1))
        sk = 0
        for n in range (Nout+1):
            for m in range (Nout+1-n):
                rout[sk] = -1 + 2*m/Nout
                sout[sk] = -1 + 2*n/Nout
                counter[n,m] = sk
                sk += 1
                
        # Build matrix to interpolate field data to equally spaced nodes
        Vout = vandermonde(Nout, rout, sout)
        interp = Vout.dot(np.linalg.inv(Vout))

        # Build triangulation of equally spaced nodes on reference triangle
        tri = np.array([], dtype=int).reshape(0,3)
        for n in range (Nout+1):
            for m in range (Nout-n):
                v1 = counter[n,m]
                v2 = counter[n,m+1]
                v3 = counter[n+1,m]
                v4 = counter[n+1,m+1]
                if v4:
                    tri = np.vstack((tri, [v1, v2, v3],[v2, v4, v3]))
                else:
                    tri = np.vstack((tri, [[v1, v2, v3]]))

        # Build triangulation for all equally spaced nodes on all elements
        TRI = np.array([], dtype=int).reshape(0,3)
        for k in range(self.mesh.number_of_elements()):
            TRI = np.vstack((TRI, tri+(k)*Npout))

        # Interpolate node coordinates and field to equally spaced nodes
        xout = interp.dot(self.x) 
        yout = interp.dot(self.y) 
        uout = interp.dot(field)

        levels = np.linspace(-1, 1, 200)
        # Render and format solution field
        plt.tricontourf(
            xout.ravel('F'), 
            yout.ravel('F'), 
            uout.ravel('F'), 
            triangles=TRI, 
            cmap='viridis',
            levels=levels
        )

    

