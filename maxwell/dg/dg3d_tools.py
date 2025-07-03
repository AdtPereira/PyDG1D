import numpy as np
import scipy

import maxwell.dg.dg1d_tools as dg1d_tools
import maxwell.dg.dg2d_tools as dg2d_tools
import maxwell.dg.mesh3d as mesh

N_FACES = 4

ALPHA_STORE = np.array([0.0000, 0.0000, 0.0000, 0.1002, 1.1332, 1.5608, 1.3413,
                        1.2577, 1.1603, 1.10153, 0.6080, 0.4523, 0.8856, 0.8717, 0.9655])

NODETOL = 1e-12

def evalwarp(p, xnodes, xout):
    """
    Compute the 1D edge warping function.

    Parameters:
    p       : int         -> Polynomial order
    xnodes  : ndarray     -> Actual nodes (usually GLL points)
    xout    : ndarray     -> Points where the warp is evaluated

    Returns:
    warp    : ndarray     -> Warping function evaluated at xout
    """

    warp = np.zeros_like(xout)

    # Equally spaced points on [-1,1]
    xeq = np.array([-1 + 2*(p + 1 - i)/p for i in range(1, p+2)])

    for i in range(p + 1):
        d = (xnodes[i] - xeq[i]) * np.ones_like(xout)
        for j in range(1, p):  # j from 2 to p in MATLAB (1 to p-1 in 0-based Python)
            if i != j:
                d *= (xout - xeq[j]) / (xeq[i] - xeq[j])
        
        if i != 0:
            d *= -1 / (xeq[i] - xeq[0])
        
        if i != p:
            d *= 1 / (xeq[i] - xeq[p])
        
        warp += d

    return warp

def evalshift(p, pval, L1, L2, L3):
    '''
        % function [dx, dy] = evalshift(p, pval, L1, L2, L3)
        % Purpose: compute two-dimensional Warp & Blend transform

        Compute the 2D Warp & Blend transformation on a triangle face.

        Parameters:
        p     : int       -> Polynomial order
        pval  : float     -> Alpha-like parameter (usually from alpopt)
        L1,L2,L3 : ndarray -> Barycentric coordinates

        Returns:
        dx, dy : ndarray -> Warp shifts in the (x, y) directions
    '''
    # 1) compute Gauss-Lobatto-Legendre node distribution
    gaussX = - dg1d_tools.jacobiGL(0, 0, p)

    # 2) compute blending function at each node for each edge
    blend1 = L2 * L3
    blend2 = L1 * L3
    blend3 = L1 * L2
    
    # 3) amount of warp for each node, for each edge
    warpfactor1 = 4 * evalwarp(p, gaussX, L3 - L2)
    warpfactor2 = 4 * evalwarp(p, gaussX, L1 - L3)
    warpfactor3 = 4 * evalwarp(p, gaussX, L2 - L1)

    # 4) combine blend & warp
    warp1 = blend1 * warpfactor1 * (1 + (pval*L1)**2)
    warp2 = blend2 * warpfactor2 * (1 + (pval*L2)**2)
    warp3 = blend3 * warpfactor3 * (1 + (pval*L3)**2)

    # 5) evaluate shift in equilateral triangle
    dx = 1.0 * warp1 + np.cos(2*np.pi/3)*warp2 + np.cos(4.0*np.pi/3)*warp3
    dy = 0.0 * warp1 + np.sin(2*np.pi/3)*warp2 + np.sin(4.0*np.pi/3)*warp3

    return dx, dy

def WarpShiftFace3D(p, pval, pval2, L1, L2, L3, L4):
    """
    Compute warp factors used in creating 3D Warp & Blend nodes on a face.

    Parameters:
    p      : int     -> Polynomial order
    pval   : ndarray -> Interpolation points on edge
    pval2  : ndarray -> Possibly unused (reserved for compatibility)
    L1,L2,L3,L4 : ndarray -> Barycentric coordinates of face points

    Returns:
    warpx, warpy : ndarray -> Warp shifts in x and y directions
    """
    dtan1, dtan2 = evalshift(p, pval, L2, L3, L4)
    warpx = dtan1
    warpy = dtan2
    return warpx, warpy

def EquiNodes3D(N):
    """
    Compute equidistributed nodes (r, s, t) on the reference tetrahedron.

    Returns:
    r, s, t : ndarray -> Coordinates in the reference tetrahedron
    """

    Np = int((N + 1)*(N + 2)*(N + 3) // 6)
    r = np.zeros(Np)
    s = np.zeros(Np)
    t = np.zeros(Np)

    sk = 0
    for n in range(1, N + 2):
        for m in range(1, N + 3 - n):
            for q in range(1, N + 4 - n - m):
                r[sk] = -1 + (q - 1) * 2 / N
                s[sk] = -1 + (m - 1) * 2 / N
                t[sk] = -1 + (n - 1) * 2 / N
                sk += 1

    return r, s, t

def set_nodes_in_equilateral_tetrahedron(N):
    '''
        Compute (x, y, z) coordinates of Warp & Blend nodes inside an equilateral tetrahedron
        for a given polynomial order N.

        % function [X,Y,Z] = Nodes3D(p)
        % Purpose:  compute Warp & Blend nodes
        % input:    p=polynomial order of interpolant
        % output:   X,Y,Z vectors of node coordinates in equilateral tetrahedron
    '''
    alpha = ALPHA_STORE[N] if N <= 14 else 1.0
    tol = 1e-10
    
    # total number of nodes in the tetrahedron
    Np = int((N+1)*(N+2)*(N+3)//6)

    # create equidistributed nodes
    # Coordenadas (r,s,t) no tetraedro de referência
    r, s, t = EquiNodes3D(N)

    # Convert (r, s, t) to barycentric coordinates
    L1 = (1 + t)/2
    L2 = (1 + s)/2
    L3 = -(1 + r + s + t)/2
    L4 = (1 + r)/2

    # set vertices of tetrahedron
    v1 = np.array([-1.0, -1/np.sqrt(3), -1/np.sqrt(6)])
    v2 = np.array([+1.0, -1/np.sqrt(3), -1/np.sqrt(6)])
    v3 = np.array([ 0.0, +2/np.sqrt(3), -1/np.sqrt(6)])
    v4 = np.array([ 0.0, +0/np.sqrt(3), +3/np.sqrt(6)])

    # orthogonal axis tangents on faces 1-4
    t1 = np.array([v2 - v1,
                   v2 - v1,
                   v3 - v2,
                   v3 - v1])
    t2 = np.array([v3 - 0.5*(v1+v2),
                   v4 - 0.5*(v1+v2),
                   v4 - 0.5*(v2+v3),
                   v4 - 0.5*(v1+v3)])
    
    # normalize tangents
    t1 /= np.linalg.norm(t1, axis=1)[:, np.newaxis]
    t2 /= np.linalg.norm(t2, axis=1)[:, np.newaxis]

    # Warp and blend for each face (accumulated in shiftXYZ)
    # form undeformed coordinates 
    # Coordenadas iniciais no tetraedro linear (sem deformação)
    XYZ = np.outer(L3, v1) + np.outer(L4, v2) + np.outer(L2, v3) + np.outer(L1, v4)
    shift = np.zeros_like(XYZ)

    for face in range(4):
        if face == 0: La, Lb, Lc, Ld = L1, L2, L3, L4
        elif face == 1: La, Lb, Lc, Ld = L2, L1, L3, L4
        elif face == 2: La, Lb, Lc, Ld = L3, L1, L4, L2
        elif face == 3: La, Lb, Lc, Ld = L4, L1, L3, L2

        # compute warp tangential to face 
        warp1, warp2 = WarpShiftFace3D(N, alpha, alpha, La, Lb, Lc, Ld)
        
        # compute volume blending
        blend = Lb * Lc * Ld

        # modify linear blend
        denom = (Lb + 0.5*La)*(Lc + 0.5*La)*(Ld + 0.5*La)

        # Aplica o blending e warp
        ids = denom > tol
        blend_mod = np.zeros_like(blend)
        blend_mod[ids] = (1 + (alpha*La[ids])**2) * blend[ids] / denom[ids]

        # compute warp & blend
        shift += np.outer(blend_mod * warp1, t1[face]) + np.outer(blend_mod * warp2, t2[face])

        # fix face warp
        corner_ids = np.where((La < tol) & ((Lb > tol) + (Lc > tol) + (Ld > tol) < 3))[0]
        shift[corner_ids, :] = np.outer(warp1[corner_ids], t1[face]) + np.outer(warp2[corner_ids], t2[face])

    XYZ += shift
    X, Y, Z = XYZ[:, 0], XYZ[:, 1], XYZ[:, 2]
    
    return X, Y, Z 

def xyz_to_rst(X, Y, Z):
    """
    % function [r,s,t] = xyztorst(x, y, z)
    % Purpose : Transfer from (x,y,z) in equilateral tetrahedron 
    %           to (r,s,t) coordinates in standard tetrahedron

    Parameters:
    X, Y, Z : ndarray
        Cartesian coordinates in the equilateral tetrahedron.

    Returns:
    r, s, t : ndarray
        Reference coordinates in the standard tetrahedron.
    """

    v1 = np.array([-1.0, -1/np.sqrt(3), -1/np.sqrt(6)])
    v2 = np.array([+1.0, -1/np.sqrt(3), -1/np.sqrt(6)])
    v3 = np.array([ 0.0, +2/np.sqrt(3), -1/np.sqrt(6)])
    v4 = np.array([ 0.0, +0/np.sqrt(3), +3/np.sqrt(6)])

    # back out right tet nodes
    # Corrige deslocamento do centro de massa para coordenadas relativas
    rhs = np.vstack([X, Y, Z]) - 0.5 * (v2 + v3 + v4 - v1)[:, np.newaxis]

    # Matriz de transformação A
    A = 0.5 * np.column_stack((v2 - v1, v3 - v1, v4 - v1))

    # Resolve sistema linear A * [r; s; t] = rhs
    RST = np.linalg.solve(A, rhs)
    r, s, t = RST[0, :], RST[1, :], RST[2, :]
    
    return r, s, t

def simplex_polynomial(a, b, c, i: int, j: int, k: int):
    '''
        % function [P] = Simplex3DP(a,b,c,i,j,k);
        % Purpose : Evaluate 3D orthonormal polynomial
        %           on simplex at (a,b,c) of order (i,j,k).
    '''
    h1 = dg1d_tools.jacobi_polynomial(a, 0, 0, i)
    h2 = dg1d_tools.jacobi_polynomial(b, 2*i + 1, 0, j)
    h3 = dg1d_tools.jacobi_polynomial(c, 2*(i + j) + 2, 0, k)

    P = np.sqrt(8.0) * h1.transpose() * h2.transpose() * (1-b)**i * h3.transpose() * (1-c)**(i+j)
    return P

def rst_to_abc(r, s, t):
    '''
        % function [a,b,c] = rst_to_abc(r,s,t)
        % Purpose : Transfer from (r,s,t) -> (a,b,c) coordinates in tetrahedron

        Parâmetros:
        -----------
        r, s, t : array_like
            Coordenadas nos elementos do tetraedro de referência.

        Retorna:
        --------
        a, b, c : ndarray
            Coordenadas transformadas.
    '''
    r = np.asarray(r)
    s = np.asarray(s)
    t = np.asarray(t)

    Np = len(r)
    a = np.zeros(Np)
    b = np.zeros(Np)
    c = np.zeros(Np)

    for n in range(Np):
        if s[n] + t[n] != 0:
            a[n] = 2 * (1 + r[n]) / (-s[n] - t[n]) - 1
        else:
            a[n] = -1

        if t[n] != 1:
            b[n] = 2 * (1 + s[n]) / (1 - t[n]) - 1
        else:
            b[n] = -1

    c = t.copy()

    return a, b, c

def vandermonde(N: int, r, s, t):
    """
    Initialize the 3D Vandermonde matrix V_{ij} = phi_j(r_i, s_i, t_i)

    Parameters
    ----------
    N : int
        Polynomial order.
    r, s, t : ndarray
        Coordinates in reference tetrahedron.

    Returns
    -------
    V3D : ndarray
        The 3D Vandermonde matrix of shape (len(r), Np), where
        Np = (N+1)(N+2)(N+3)/6
    """
    Np = int((N+1)*(N+2)*(N+3)//6)
    V3D = np.zeros((len(r), Np))

    # Convert (r, s, t) to barycentric coordinates (a, b, c)
    a, b, c = rst_to_abc(r, s, t)

    # Build the Vandermonde matrix
    sk = 0
    for i in range(N + 1):
        for j in range(N + 1 - i):
            for k in range(N + 1 - i - j):
                V3D[:, sk] = simplex_polynomial(a, b, c, i, j, k)
                sk += 1

    return V3D

def massMatrix(N, r, s, t):
    vander = vandermonde(N, r, s, t)
    mass = np.linalg.inv(vander.dot(vander.transpose()))
    return mass

def derivateMatrix(N, r, s, t, V):
    """
    % function [Dr,Ds,Dt] = Dmatrices3D(N,r,s,t,V)
    % Purpose : Initialize the (r,s,t) differentiation matrices
    %           on the simplex, evaluated at (r,s,t) at order N

    Inicializa as matrizes de diferenciação (Dr, Ds, Dt) no tetraedro de referência.

    Parâmetros:
    -----------
    N : int
        Ordem do polinômio.
    r, s, t : ndarray
        Coordenadas nos pontos nodais.
    V : ndarray
        Matriz de Vandermonde 3D avaliada em (r,s,t).

    Retorna:
    --------
    Dr, Ds, Dt : ndarray
        Matrizes de derivadas nas direções r, s, t.
    """
    Vr, Vs, Vt = gradVandermonde(N, r, s, t)

    # Matrizes diferenciais: Dr = Vr @ inv(V)
    Vinv = np.linalg.inv(V)
    Dr = Vr @ Vinv
    Ds = Vs @ Vinv
    Dt = Vt @ Vinv

    return Dr, Ds, Dt

def gradVandermonde(N: int, r, s, t):
    """
    % function [V3Dr,V3Ds,V3Dt] = GradVandermonde3D(N,r,s,t)
    % Purpose : Initialize the gradient of the modal basis (i,j,k) 
    %           at (r,s,t) at order N
    Inicializa as derivadas da base modal (i,j,k) nos pontos (r,s,t) até a ordem N.

    Parâmetros:
    -----------
    N : int
        Ordem do polinômio.
    r, s, t : ndarray
        Coordenadas no tetraedro de referência.

    Retorna:
    --------
    V3Dr, V3Ds, V3Dt : ndarray
        Matrizes com as derivadas das funções de base modais nas direções r, s e t.
    """

    Np = int((N+1)*(N+2)*(N+3)//6)
    V3Dr = np.zeros((len(r), Np))
    V3Ds = np.zeros((len(r), Np))
    V3Dt = np.zeros((len(r), Np))

    # find tensor-product coordinates
    a, b, c = rst_to_abc(r,s,t)

    # Initialize matrices
    sk = 0
    for i in range(0, N + 1, 1):
        for j in range(0, N - i + 1, 1):
            for k in range(0, N -i - j + 1, 1):
                V3Dr[:, sk], V3Ds[:, sk], V3Dt[:, sk] = gradSimplexP(a, b, c, i, j, k)
                sk += 1
    
    return V3Dr, V3Ds, V3Dt

def gradSimplexP(a, b, c, id: int, jd: int, kd: int):
    """
    % function [V3Dr, V3Ds, V3Dt] = GradSimplex3DP(a,b,c,id,jd,kd)
    % Purpose: Return the derivatives of the modal basis (id,jd,kd)
    %          on the 3D simplex at (a,b,c)

    Calcula as derivadas das funções de base modais no tetraedro de referência.

    Parâmetros:
    -----------
    a, b, c : ndarray
        Coordenadas transformadas (a,b,c) no tetraedro.
    id, jd, kd : int
        Índices do polinômio de Jacobi na base modal.

    Retorna:
    --------
    V3Dr, V3Ds, V3Dt : ndarray
        Derivadas das funções de base modais nas direções r, s e t.
    """

    # Polinômios de Jacobi e suas derivadas
    fa  = dg1d_tools.jacobi_polynomial     (a, 0, 0, id)
    dfa = dg1d_tools.jacobi_polynomial_grad(a, 0, 0, id)
    gb  = dg1d_tools.jacobi_polynomial     (b, 2.0*id + 1.0, 0, jd)
    dgb = dg1d_tools.jacobi_polynomial_grad(b, 2.0*id + 1.0, 0, jd)
    hc  = dg1d_tools.jacobi_polynomial     (c, 2.0*(id+jd) + 2.0, 0, kd)
    dhc = dg1d_tools.jacobi_polynomial_grad(c, 2.0*(id+jd) + 2.0, 0, kd)

    # r-derivative
    V3Dr = dfa * gb * hc
    if id > 0:
        V3Dr *= (0.5 * (1 - b)) ** (id - 1)
    if id + jd > 0:
        V3Dr *= (0.5 * (1 - c)) ** (id + jd - 1)

    # s-derivative
    V3Ds = 0.5 * (1 + a) * V3Dr
    tmp = dgb * (0.5 * (1 - b)) ** id
    if id > 0:
        tmp += -0.5 * id * gb * (0.5 * (1 - b)) ** (id - 1)
    if id + jd > 0:
        tmp *= (0.5 * (1 - c)) ** (id + jd - 1)
    tmp = fa * tmp * hc
    V3Ds += tmp

    # t-derivative
    V3Dt = 0.5 * (1 + a) * V3Dr + 0.5 * (1 + b) * tmp
    tmp = dhc * (0.5 * (1 - c)) ** (id + jd)
    if id + jd > 0:
        tmp -= 0.5 * (id + jd) * hc * (0.5 * (1 - c)) ** (id + jd - 1)
    tmp = fa * gb * tmp * (0.5 * (1 - b)) ** id
    V3Dt += tmp

    # Normalização
    factor = 2 ** (2 * id + jd + 1.5)
    V3Dr *= factor
    V3Ds *= factor
    V3Dt *= factor

    return V3Dr, V3Ds, V3Dt

def map_reference_to_physical_nodes(msh, r_ref, s_ref, t_ref):
    """
    Mapeia um conjunto de nós (r, s, t) do tetraedro de referência para as
    coordenadas físicas (x, y, z) de todos os elementos da malha.

    Parâmetros
    ----------
    msh : Mesh3D
        O objeto de malha contendo a conectividade e as coordenadas dos vértices.
    r_ref, s_ref, t_ref : ndarray
        Arrays 1D com as coordenadas dos nós no elemento de referência.

    Retorna
    -------
    x, y, z : ndarray
        Arrays com as coordenadas físicas dos nós, com shape (len(r_ref), K).
    """
    # Índices dos vértices para cada elemento.
    va = msh.EToV[:, 0]
    vb = msh.EToV[:, 1]
    vc = msh.EToV[:, 2]
    vd = msh.EToV[:, 3]

    # Garante que os arrays de referência tenham a forma (Np, 1) para broadcasting.
    r = r_ref.reshape(-1, 1)
    s = s_ref.reshape(-1, 1)
    t = t_ref.reshape(-1, 1)

    # Fórmula de interpolação linear baseada nos vértices de cada tetraedro.
    x = 0.5 * (-(1 + r + s + t) * msh.vx[va] + (1 + r) * msh.vx[vb] + (1 + s) * msh.vx[vc] + (1 + t) * msh.vx[vd])
    y = 0.5 * (-(1 + r + s + t) * msh.vy[va] + (1 + r) * msh.vy[vb] + (1 + s) * msh.vy[vc] + (1 + t) * msh.vy[vd])
    z = 0.5 * (-(1 + r + s + t) * msh.vz[va] + (1 + r) * msh.vz[vb] + (1 + s) * msh.vz[vc] + (1 + t) * msh.vz[vd])

    return x, y, z

def nodesCoordinates(N, msh: mesh.Mesh3D):
    """
    Constrói as coordenadas físicas (x, y, z) dos nós de solução DG
    para uma dada ordem N.
    """
    # Gera os nós de referência específicos para a ordem N.
    r, s, t = xyz_to_rst(*set_nodes_in_equilateral_tetrahedron(N))
    
    # Chama a função de mapeamento genérica.
    return map_reference_to_physical_nodes(msh, r, s, t)

def buildFMask(N):
    r, s, t  = xyz_to_rst(*set_nodes_in_equilateral_tetrahedron(N))
    fmask1 = np.where(np.abs(t+1) < NODETOL)[0]
    fmask2 = np.where(np.abs(s+1) < NODETOL)[0]
    fmask3 = np.where(np.abs(r+s+t+1) < NODETOL)[0]
    fmask4 = np.where(np.abs(r+1) < NODETOL)[0]
    Fmask  = np.array([fmask1, fmask2, fmask3, fmask4]).transpose()

    return Fmask, fmask1, fmask2, fmask3, fmask4

def lift(N):
    """
    % function [LIFT] = Lift3D(N, r, s, t)
    % Purpose : Compute 3D surface to volume lift operator used 
    %           in DG formulation

    Calcula o operador de elevação (LIFT) para termos de superfície no método DG 3D.

    Parâmetros:
        N : int
            Ordem do polinômio.
        R, S, T : ndarray (Np,)
            Coordenadas dos nós no tetraedro de referência.
        V : ndarray (Np, Np)
            Matriz de Vandermonde 3D no tetraedro.
        Fmask : ndarray (Nfp, Nfaces)
            Máscara com índices dos nós em cada face.
        Nfaces : int
            Número de faces do tetraedro (default: 4).

    Retorna:
        LIFT : ndarray (Np, Nfaces*Nfp)
            Operador de elevação.
    """
    Nfp = int((N+1)*(N+2)//2)
    Np = int((N+1)*(N+2)*(N+3)//6)
    Emat = np.zeros((Np, N_FACES*Nfp))

    r, s, t = xyz_to_rst(*set_nodes_in_equilateral_tetrahedron(N))
    Fmask, _, _, _, _ = buildFMask(N)
    
    for face in range(N_FACES):
        ids = Fmask[:, face]
        if face == 0:
            faceR, faceS = r[ids], s[ids]
        if face == 1:
            faceR, faceS = r[ids], t[ids]
        if face == 2:
            faceR, faceS = s[ids], t[ids]
        if face == 3:
            faceR, faceS = s[ids], t[ids]

        VFace = dg2d_tools.vandermonde(N, faceR, faceS)
        massFace = np.linalg.inv(VFace @ VFace.T)

        idc = slice(face * Nfp, (face + 1) * Nfp)
        Emat[ids, idc] = Emat[ids, idc] + massFace

    # LIFT = inv(MassMatrix) @ Emat, com MassMatrix ≈ V.T @ V
    V = vandermonde(N, r, s, t)
    return V @ (V.T @ Emat)

def geometricFactors(x, y, z, Dr, Ds, Dt):
    """
    % function [rx,sx,tx,ry,sy,ty,rz,sz,tz,J] = ...
    % GeometricFactors3D(x,y,z,Dr,Ds,Dt)
    % Purpose : Compute the metric elements for the local mappings of the elements

    Calcula os fatores geométricos para o mapeamento local dos elementos tetraédricos.

    Parâmetros:
        x, y, z : ndarray (Np, K)
            Coordenadas físicas dos nós em cada elemento.
        Dr, Ds, Dt : ndarray (Np, Np)
            Matrizes de diferenciação nas direções r, s e t.

    Retorna:
        rx, sx, tx,
        ry, sy, ty,
        rz, sz, tz : ndarray (Np, K)
            Derivadas das coordenadas físicas em relação às coordenadas referenciais.
        J : ndarray (Np, K)
            Jacobiano do mapeamento.
    """
    # calculate geometric factors
    xr = Dr @ x
    xs = Ds @ x
    xt = Dt @ x

    yr = Dr @ y
    ys = Ds @ y
    yt = Dt @ y

    zr = Dr @ z
    zs = Ds @ z
    zt = Dt @ z

    # Jacobian determinant
    J = xr * (ys * zt - zs * yt) \
        - yr * (xs * zt - zs * xt) \
        + zr * (xs * yt - ys * xt)

    # Inverse mapping factors
    rx =  (ys * zt - zs * yt) / J
    ry = -(xs * zt - zs * xt) / J
    rz =  (xs * yt - ys * xt) / J

    sx = -(yr * zt - zr * yt) / J
    sy =  (xr * zt - zr * xt) / J
    sz = -(xr * yt - yr * xt) / J

    tx =  (yr * zs - zr * ys) / J
    ty = -(xr * zs - zr * xs) / J
    tz =  (xr * ys - yr * xs) / J

    return rx, sx, tx, ry, sy, ty, rz, sz, tz, J

def normals(x, y, z, Dr, Ds, Dt, N):
    """
    % function [nx, ny, nz, sJ] = Normals3D()
    % Purpose : Compute outward pointing normals 
    %           at elements faces % as well as surface Jacobians.

    Calcula os vetores normais externos nos nós das faces e os Jacobianos de superfície sJ.

    Parâmetros:
        x, y, z : ndarray (Np, K)
            Coordenadas físicas dos nós.
        Dr, Ds, Dt : ndarray (Np, Np)
            Matrizes diferenciais no referencial.
        Fmask : ndarray (Nfp, 4)
            Índices dos nós de cada face do elemento.
        Nfp : int
            Número de nós por face.

    Retorna:
        nx, ny, nz : ndarray (4*Nfp, K)
            Componentes normalizadas das normais externas por face.
        sJ : ndarray (4*Nfp, K)
            Jacobianos de superfície multiplicados por |n|.
    """
    Nfp = int((N+1)*(N+2)//2)
    Fmask = buildFMask(N)[0]
    rx, sx, tx, ry, sy, ty, rz, sz, tz, J = geometricFactors(x, y, z, Dr, Ds, Dt)

    # Interpola os fatores geométricos para os nós de face
    ids = Fmask.flatten()
    frx = rx[ids, :]
    fsx = sx[ids, :]
    ftx = tx[ids, :]
    fry = ry[ids, :]
    fsy = sy[ids, :]
    fty = ty[ids, :]
    frz = rz[ids, :]
    fsz = sz[ids, :]
    ftz = tz[ids, :]

    # build normals
    K = x.shape[1]
    nx = np.zeros((4 * Nfp, K))
    ny = np.zeros((4 * Nfp, K))
    nz = np.zeros((4 * Nfp, K))

    fid1 = np.arange(0, Nfp)
    fid2 = np.arange(Nfp, 2 * Nfp)
    fid3 = np.arange(2 * Nfp, 3 * Nfp)
    fid4 = np.arange(3 * Nfp, 4 * Nfp)

    # Face 1
    nx[fid1, :] = -ftx[fid1, :]
    ny[fid1, :] = -fty[fid1, :]
    nz[fid1, :] = -ftz[fid1, :]

    # Face 2
    nx[fid2, :] = -fsx[fid2, :]
    ny[fid2, :] = -fsy[fid2, :]
    nz[fid2, :] = -fsz[fid2, :]

    # Face 3
    nx[fid3, :] = frx[fid3, :] + fsx[fid3, :] + ftx[fid3, :]
    ny[fid3, :] = fry[fid3, :] + fsy[fid3, :] + fty[fid3, :]
    nz[fid3, :] = frz[fid3, :] + fsz[fid3, :] + ftz[fid3, :]

    # Face 4
    nx[fid4, :] = -frx[fid4, :]
    ny[fid4, :] = -fry[fid4, :]
    nz[fid4, :] = -frz[fid4, :]

    # Normalização
    sJ = np.sqrt(nx**2 + ny**2 + nz**2)
    nx /= sJ
    ny /= sJ
    nz /= sJ

    sJ *= J[ids, :]

    return nx, ny, nz, sJ

def Grad3D(U, Dr, Ds, Dt, rx, sx, tx, ry, sy, ty, rz, sz, tz):
    """
    % function [dUdx, dUdy, dUdz] = GradH3D(U)
    % purpose: Compute local elemental physical spatial 
    %           derivatives of U in 3D.

    Parameters
    ----------
    U : ndarray
        Nodal values of the scalar field (shape: [Np, K]).
    Dr, Ds, Dt : ndarray
        Differentiation matrices in r, s, t directions.
    rx, sx, tx : ndarray
        Metric terms for x-derivatives (shape: [Np, K]).
    ry, sy, ty : ndarray
        Metric terms for y-derivatives.
    rz, sz, tz : ndarray
        Metric terms for z-derivatives.

    Returns
    -------
    dUdx, dUdy, dUdz : ndarray
        Derivatives of U in physical space.
    """
    # Reference derivatives
    dUdr = Dr @ U
    dUds = Ds @ U
    dUdt = Dt @ U

    # Physical derivatives using the chain rule
    dUdx = rx * dUdr + sx * dUds + tx * dUdt
    dUdy = ry * dUdr + sy * dUds + ty * dUdt
    dUdz = rz * dUdr + sz * dUds + tz * dUdt

    return dUdx, dUdy, dUdz

def Div3D(Ux, Uy, Uz, Dr, Ds, Dt, rx, sx, tx, ry, sy, ty, rz, sz, tz):
    """
    % function [divU] = DivH3D(Ux, Uy, Uz)
    % Purpose: Compute local elemental physical spatial divergence of (Ux, Uy, Uz) in 3D.

    Parameters
    ----------
    Ux, Uy, Uz : ndarray
        Components of the vector field (shape: [Np, K]).
    Dr, Ds, Dt : ndarray
        Differentiation matrices in r, s, t directions.
    rx, sx, tx : ndarray
        Metric terms for x-derivatives.
    ry, sy, ty : ndarray
        Metric terms for y-derivatives.
    rz, sz, tz : ndarray
        Metric terms for z-derivatives.

    Returns
    -------
    divU : ndarray
        Divergence of the vector field (shape: [Np, K]).
    """
    # compute local derivatives of Ux on reference tetrahedron
    ddr = Dr @ Ux
    dds = Ds @ Ux
    ddt = Dt @ Ux
    divU = rx * ddr + sx * dds + tx * ddt

    # compute local derivatives of Uy on reference tetrahedron
    ddr = Dr @ Uy
    dds = Ds @ Uy
    ddt = Dt @ Uy
    divU += ry * ddr + sy * dds + ty * ddt

    # compute local derivatives of Uz on reference tetrahedron
    ddr = Dr @ Uz
    dds = Ds @ Uz
    ddt = Dt @ Uz
    divU += rz * ddr + sz * dds + tz * ddt

    return divU

def Curl3D(Ux, Uy, Uz, Dr, Ds, Dt, rx, sx, tx, ry, sy, ty, rz, sz, tz):
    """
    % function [curlx, curly, curlz] = Curl3D(Ux, Uy, Uz)
    % purpose: compute local elemental physical spatial curl of (Ux,Uy,Uz) in 3D.

    Parameters
    ----------
    Ux, Uy, Uz : ndarray
        Components of the vector field (shape: [Np, K]).
    Dr, Ds, Dt : ndarray
        Derivative matrices with respect to (r, s, t).
    rx, sx, tx : ndarray
        Metric terms for ∂/∂x.
    ry, sy, ty : ndarray
        Metric terms for ∂/∂y.
    rz, sz, tz : ndarray
        Metric terms for ∂/∂z.

    Returns
    -------
    curlx, curly, curlz : ndarray
        Components of the curl of the field.
    """

    # compute local derivatives of Ux on reference tetrahedron
    ddr = Dr @ Ux
    dds = Ds @ Ux
    ddt = Dt @ Ux

    # increment curl components
    curly = rz * ddr + sz * dds + tz * ddt
    curlz = - (ry * ddr + sy * dds + ty * ddt)

    # compute local derivatives of Uy on reference tetrahedron
    ddr = Dr @ Uy
    dds = Ds @ Uy
    ddt = Dt @ Uy

    # increment curl components
    curlx = - (rz * ddr + sz * dds + tz * ddt)
    curlz += rx * ddr + sx * dds + tx * ddt

    # compute local derivatives of Uz on reference tetrahedron
    ddr = Dr @ Uz
    dds = Ds @ Uz
    ddt = Dt @ Uz

    # increment curl components
    curlx += ry * ddr + sy * dds + ty * ddt
    curly -= rx * ddr + sx * dds + tx * ddt

    return curlx, curly, curlz

