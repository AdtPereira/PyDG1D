import numpy as np
from scipy.linalg import eigh
from scipy.special import gamma
import matplotlib.pyplot as plt
from matplotlib.tri import Triangulation
import pandas as pd
from scipy.sparse import lil_matrix, eye
import scipy.special
from math import ceil, pi
from mpl_toolkits.mplot3d import Axes3D
from meshpy.tet import MeshInfo, build
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
from mpl_toolkits.mplot3d.art3d import Poly3DCollection, Line3DCollection
import vtk
from scipy.interpolate import griddata
import matplotlib.patches as patches
import gmsh
import os

# Limpar terminal
os.system('cls' if os.name == 'nt' else 'clear')

def JacobiP(x, alpha, beta, N):
    """
    Avalia o polinômio de Jacobi ortonormal P_N^{(alpha,beta)} em x.
    
    Parâmetros
    ----------
    x : array-like
        Pontos em que avaliar o polinômio.
    alpha, beta : floats
        Parâmetros do polinômio (devem ser > -1 e alpha+beta != -1).
    N : int
        Grau do polinômio.
    
    Retorna
    -------
    P : ndarray
        Vetor com P_N^{(alpha,beta)}(x) para cada x.
    """
    x = np.asarray(x, dtype=float)
    # Garantir array 1D
    xp = x.ravel()
    M = xp.size

    # Matriz para P_0 ... P_N
    PL = np.zeros((N+1, M))

    # P_0
    gamma0 = 2**(alpha+beta+1) / (alpha+beta+1) * gamma(alpha+1) * gamma(beta+1) / gamma(alpha+beta+1)
    PL[0, :] = 1.0 / np.sqrt(gamma0)
    if N == 0:
        return PL[0, :].reshape(x.shape)

    # P_1
    gamma1 = (alpha+1)*(beta+1)/(alpha+beta+3) * gamma0
    PL[1, :] = (((alpha+beta+2) * xp / 2) + (alpha-beta)/2) / np.sqrt(gamma1)
    if N == 1:
        return PL[1, :].reshape(x.shape)

    # Recorrência
    aold = 2 / (2 + alpha + beta) * np.sqrt((alpha+1)*(beta+1)/(alpha+beta+3))
    for i in range(1, N):
        h1 = 2*i + alpha + beta
        anew = (2 / (h1 + 2) *
                np.sqrt((i+1)*(i+1+alpha+beta)*(i+1+alpha)*(i+1+beta)/
                        ((h1+1)*(h1+3))))
        bnew = - (alpha**2 - beta**2) / (h1 * (h1+2))
        PL[i+1, :] = ( (xp - bnew) * PL[i, :] - aold * PL[i-1, :] ) / anew
        aold = anew

    # Extrai P_N e reformata para a forma original de x
    PN = PL[N, :].reshape(x.shape)
    return PN

def Simplex3DP(a, b, c, i, j, k):
    """
    Evaluate 3D orthonormal polynomial on a simplex at (a, b, c) of order (i, j, k).

    Parameters
    ----------
    a, b, c : array_like or float
        Coordinates in the reference simplex (each can be a NumPy array or scalar).
    i, j, k : int
        Polynomial orders in the a, b, and c directions.

    Returns
    -------
    P : ndarray or float
        The value(s) of the 3D orthonormal polynomial at (a, b, c).
    """
    # Evaluate 1D Jacobi polynomials
    h1 = JacobiP(a, 0, 0, i)
    h2 = JacobiP(b, 2*i + 1, 0, j)
    h3 = JacobiP(c, 2*(i + j) + 2, 0, k)

    # Build the 3D simplex polynomial
    P = 2 * np.sqrt(2) * h1 * h2 * (1 - b)**i * h3 * (1 - c)**(i + j)
    return P

def rsttoabc(r, s, t):

    """
    Transforma coordenadas (r, s, t) em coordenadas (a, b, c) no triângulo.

    Parâmetros
    ----------
    r, s, t : array-like, shape (Np,)
        Coordenadas de entrada.

    Retorna
    -------
    a, b, c : ndarray, shape (Np,)
        Coordenadas transformadas.
    """
    r = np.asarray(r)
    s = np.asarray(s)
    t = np.asarray(t)

    # Inicializa com valor default
    a = np.full_like(r, -1, dtype=float)
    b = np.full_like(r, -1, dtype=float)

    # Para os índices onde s + t != 0
    mask_a = (s + t) != 0
    a[mask_a] = 2 * (1 + r[mask_a]) / (-(s[mask_a] + t[mask_a])) - 1

    # Para os índices onde t != 1
    mask_b = t != 1
    b[mask_b] = 2 * (1 + s[mask_b]) / (1 - t[mask_b]) - 1

    # c é simplesmente t
    c = t.copy()

    return a, b, c

def JacobiGQ(alpha: float, beta: float, N: int):
    '''
    Função para computar os pontos e pesos de gauss para quadratura gaussiana
    no intervalo de [-1, 1]
    '''

    n = np.arange(N)


    # 1) Diagonal: a_n
    if abs(alpha + beta) < 1e-14:
        a = np.zeros(N)  # caso Legendre puro
    else:
        a = (beta**2 - alpha**2) / ((2*n + alpha + beta) * (2*n + alpha + beta + 2))

    # 2) Sub/super-diagonal: b_n
    b = np.zeros(N-1)
    for i in range(N-1):
        num = (i+1)*(i+alpha+1)*(i+beta+1)*(i+alpha+beta+1)
        den = (2*i+alpha+beta+1)*(2*i+alpha+beta+3)
        b[i] = 2/(2*i+alpha+beta+2) * np.sqrt(num/den)

    # 3) Monta a matriz de Jacobi tridiagonal
    J = np.diag(a) + np.diag(b,1) + np.diag(b,-1)

    # 4) Autovalores (x) e autovetores (V) de J
    x, V = eigh(J)

    # 5) Pesos
    prefac = 2**(alpha+beta+1) * gamma(alpha+1) * gamma(beta+1) / gamma(alpha+beta+2)
    w = prefac * (V[0,:]**2)

    return x, w

def JacobiGL(alpha, beta, n_order):
    """
    Compute the order N Gauss Lobatto quadrature points, x, associated
    with the Jacobi polynomial.

    >>> jacobi_gauss_lobatto(0.0, 0.0, 1)
    array([-1.,  1.])

    >>> jacobi_gauss_lobatto(0,0,3)
    array([-1.       , -0.4472136,  0.4472136,  1.       ])

    >>> jacobi_gauss_lobatto(0,0,4)
    array([-1.        , -0.65465367,  0.        ,  0.65465367,  1.        ])

    """
    if n_order == 0:
        return np.array([0.0])
    if n_order == 1:
        return np.array([-1.0, 1.0])
    if n_order > 1:
        x, w = scipy.special.roots_jacobi(n_order-1, alpha+1, beta+1)
        r_gl = np.concatenate(([-1.0], x, [1.0]))
        return r_gl

    raise ValueError('n_order must be positive.')

def GradJacobiP(r, alpha, beta, N):
    """
    Avalia a derivada do polinômio de Jacobi de tipo (alpha, beta) no(s) ponto(s) r,
    para ordem N.

    Parâmetros
    ----------
    r : np.ndarray
        Ponto(s) de avaliação, array unidimensional.
    alpha : float
        Parâmetro alpha > -1.
    beta : float
        Parâmetro beta > -1.
    N : int
        Grau do polinômio.

    Retorna
    -------
    dP : np.ndarray
        Valores da derivada P_N'^(alpha,beta)(r), mesmo formato de `r`.
    """
    r = np.asarray(r)
    dP = np.zeros_like(r, dtype=float)

    if N == 0:
        # Para N=0, a derivada é zero em todos os pontos
        return dP

    # A relação: d/dx P_N^(alpha,beta)(x) = sqrt(N*(N+alpha+beta+1)) * P_{N-1}^{(alpha+1,beta+1)}(x)
    dP = np.sqrt(N * (N + alpha + beta + 1)) * JacobiP(r[:], alpha+1, beta+1, N-1)

    return dP

def evalwarp(p, xnodes, xout):
    """
    Computa a função de warping unidimensional usada no Warp & Blend.

    Parâmetros
    ----------
    p : int
        Ordem do polinômio (de 1 em diante).
    xnodes : array-like, shape (p+1,)
        Nós de Gauss–Lobatto–Legendre para uma aresta.
    xout : array-like
        Pontos de saída onde calcular o warping.

    Retorna
    -------
    warp : ndarray, shape same as xout
        Valores do warping em cada ponto de xout.
    """
    warp = np.zeros_like(xout, dtype=float)

    # Coordenadas de um triângulo equidistante: xeq[0]=+1, xeq[p]=-1
    xeq = np.array([-1 + 2*(p + 1 - i)/p for i in range(1, p+2)])

    for i in range(p + 1):
        # Termo inicial: diferença entre nó e posição equidistante
        d = xnodes[i] - xeq[i]

        # Produto sobre j = 2..p (em MATLAB) corresponde a j_idx = 1..p-1 em Python
        for j_idx in range(1, p):
            if j_idx != i:
                d = d * (xout - xeq[j_idx]) / (xeq[i] - xeq[j_idx])

        # Ajuste de borda esquerda (i != 1 em MATLAB → i_idx != 0)
        if i != 0:
            d = -d / (xeq[i] - xeq[0])

        # Ajuste de borda direita (i != p+1 em MATLAB → i_idx != p)
        if i != p:
            d = d / (xeq[i] - xeq[p])

        warp += d

    return warp

def evalshift(p, pval, L1, L2, L3):
    """
    Computa o deslocamento (dx, dy) do Warp & Blend para um triângulo.

    Parâmetros
    ----------
    p : array-like, shape (Np,)
        Coordenada 1D de cada nó (ex.: r ou s).
    pval : array-like, shape (Np,)
        Valor auxiliar associado a cada nó (normalmente coord. de mistura).
    L1, L2, L3 : array-like, shape (Np,)
        Coordenadas barycêntricas de cada nó.

    Retorna
    -------
    dx, dy : ndarray, shape (Np,)
        Deslocamentos resultantes na direção x e y para cada nó.
    """
    

    # 1) distribuição Gauss–Lobatto–Legendre
    gaussX = - JacobiGL(0, 0, p)

    # 2) funções de blend para cada aresta
    blend1 = L2 * L3
    blend2 = L1 * L3
    blend3 = L1 * L2

    # 3) fator de warp para cada aresta
    warpfactor1 = 4 * evalwarp(p, gaussX, L3 - L2)
    warpfactor2 = 4 * evalwarp(p, gaussX, L1 - L3)
    warpfactor3 = 4 * evalwarp(p, gaussX, L2 - L1)

    # 4) combinação de blend e warp
    warp1 = blend1 * warpfactor1 * (1 + (pval * L1) ** 2)
    warp2 = blend2 * warpfactor2 * (1 + (pval * L2) ** 2)
    warp3 = blend3 * warpfactor3 * (1 + (pval * L3) ** 2)

    # 5) projeção para a transformação no triângulo equilátero
    dx = (
        warp1
        + np.cos(2 * np.pi / 3) * warp2
        + np.cos(4 * np.pi / 3) * warp3
    )
    dy = (
        0 * warp1
        + np.sin(2 * np.pi / 3) * warp2
        + np.sin(4 * np.pi / 3) * warp3
    )

    return dx, dy

def WarpShiftFace3D(p, pval, pval2, L1, L2, L3, L4):
    """
    Computa o fator de warp usado na criação de nós 3D Warp & Blend.

    Parâmetros
    ----------
    p : array-like, shape (Np, 3)
        Coordenadas dos pontos.
    pval : array-like
        Valores de referência para o warp.
    pval2 : array-like
        Segundo conjunto de valores de referência (não usado nesta função mas mantido na assinatura).
    L1, L2, L3, L4 : float
        Parâmetros de escala para a transformação.

    Retorna
    -------
    warpx, warpy : ndarray
        Componentes do deslocamento de warp em x e y.
    """

    # Chama evalshift para obter dtan1, dtan2
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

    Np = (N + 1)*(N + 2)*(N + 3) // 6
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

def Nodes3D(p):
    """
    Compute Warp & Blend nodes for an equilateral tetrahedron of order p.

    Parameters
    ----------
    p : int
        Polynomial order of the interpolant.

    Returns
    -------
    X, Y, Z : ndarray of shape (N,)
        Coordinates of the N = (p+1)*(p+2)*(p+3)/6 nodes.
    """
    # Optimized blending parameters for p up to 15
    alphastore = np.array([
        0.0, 0.0, 0.0, 0.1002, 1.1332, 1.5608,
        1.3413, 1.2577, 1.1603, 1.10153, 0.6080,
        0.4523, 0.8856, 0.8717, 0.9655
    ])
    alpha = float(alphastore[p]) if p <= 14 else 1.0

    # Number of nodes and tolerance
    N = int((p+1)*(p+2)*(p+3)//6)
    tol = 1e-10

    # Equidistributed barycentric coordinates
    r, s, t = EquiNodes3D(p)  # each of shape (N,)

    # Barycentric combinations
    L1 = (1 + t) / 2
    L2 = (1 + s) / 2
    L3 = -(1 + r + s + t) / 2
    L4 = (1 + r) / 2

    # Vertex coordinates of the equilateral tetrahedron
    v1 = np.array([-1.0, -1/np.sqrt(3), -1/np.sqrt(6)])
    v2 = np.array([ 1.0, -1/np.sqrt(3), -1/np.sqrt(6)])
    v3 = np.array([ 0.0,  2/np.sqrt(3), -1/np.sqrt(6)])
    v4 = np.array([ 0.0,  0.0,         3/np.sqrt(6)])

    # Compute orthogonal tangents on each face
    t1 = np.zeros((4, 3))
    t2 = np.zeros((4, 3))
    # Face 1 (v1,v2,v3)
    t1[0] = v2 - v1
    t1[2] = v3 - v2
    t1[3] = v3 - v1
    # Face 2 (v1,v2,v4)
    t1[1] = v2 - v1
    t2[0] = v3 - 0.5 * (v1 + v2)
    t2[1] = v4 - 0.5 * (v1 + v2)
    t2[2] = v4 - 0.5 * (v2 + v3)
    t2[3] = v4 - 0.5 * (v1 + v3)

    # Normalize tangents
    for n in range(4):
        t1[n] /= np.linalg.norm(t1[n])
        t2[n] /= np.linalg.norm(t2[n])

    # Undeformed coordinates in 3D
    XYZ = (L3[:, None] * v1 + L4[:, None] * v2
           + L2[:, None] * v3 + L1[:, None] * v4)

    # Initialize shift array
    shift = np.zeros_like(XYZ)

    # Warp & blend on each face
    for face in range(4):
        # Select barycentric labels for this face
        if face == 0:
            La, Lb, Lc, Ld = L1, L2, L3, L4
        elif face == 1:
            La, Lb, Lc, Ld = L2, L1, L3, L4
        elif face == 2:
            La, Lb, Lc, Ld = L3, L1, L4, L2
        else:
            La, Lb, Lc, Ld = L4, L1, L3, L2

        # Compute warp displacements on this face
        warp1, warp2 = WarpShiftFace3D(p, alpha, alpha, La, Lb, Lc, Ld)

        # Volume blending factor
        blend = Lb * Lc * Ld
        denom = (Lb + 0.5 * La) * (Lc + 0.5 * La) * (Ld + 0.5 * La)
        idx = denom > tol
        blend[idx] = (1 + (alpha * La[idx])**2) * blend[idx] / denom[idx]

        # Accumulate shift due to warp & blend
        shift += np.outer(blend * warp1, t1[face])
        shift += np.outer(blend * warp2, t2[face])

        # Correct pure face warp on the triangular face
        face_idx = (La < tol) & (((Lb > tol).astype(int)
                    + (Lc > tol).astype(int)
                    + (Ld > tol).astype(int)) < 3)
        shift[face_idx] = np.outer(warp1[face_idx], t1[face])
        shift[face_idx] += np.outer(warp2[face_idx], t2[face])

    # Apply the shift deformation
    XYZ += shift

    # Return separate coordinate arrays
    X = XYZ[:, 0]
    Y = XYZ[:, 1]
    Z = XYZ[:, 2]
    return X, Y, Z

def xyztorst(X, Y, Z):
    """
    Transfer from (x, y, z) in equilateral tetrahedron
    to (r, s, t) coordinates in standard tetrahedron.

    Parameters
    ----------
    X, Y, Z : array_like, shape (N,)
        Coordinates in the equilateral tetrahedron.

    Returns
    -------
    r, s, t : ndarray, shape (N,)
        Corresponding coordinates in the standard tetrahedron.
    """

    # Define equilateral tetrahedron vertex vectors
    v1 = np.array([-1.0, -1.0/np.sqrt(3), -1.0/np.sqrt(6)])
    v2 = np.array([ 1.0, -1.0/np.sqrt(3), -1.0/np.sqrt(6)])
    v3 = np.array([ 0.0,  2.0/np.sqrt(3), -1.0/np.sqrt(6)])
    v4 = np.array([ 0.0,  0.0,           3.0/np.sqrt(6)])

    # Stack input coordinates into a 3xN array
    XYZ = np.vstack((X, Y, Z))  # shape (3, N)

    # Compute the right-hand side: XYZ - midpoint adjustment
    shift = 0.5 * (v2 + v3 + v4 - v1)[:, None]  # shape (3, 1)
    rhs = XYZ - shift  # shape (3, N)

    # Build the transformation matrix A (3x3)
    A = np.column_stack((0.5 * (v2 - v1),
                         0.5 * (v3 - v1),
                         0.5 * (v4 - v1)))  # shape (3, 3)

    # Solve for RST in A * RST = rhs for each column
    RST = np.linalg.solve(A, rhs)  # shape (3, N)

    # Extract r, s, t and return as column vectors
    r = RST[0, :]
    s = RST[1, :]
    t = RST[2, :]

    return r, s, t

def Simplex2DP(a, b, i, j):

    """
    Calcula o polinômio ortonormal 2D em um simplex no ponto (a,b) para ordem (i, j)

    Parâmetros
    ----------
    a : array_like or float
        First coordinate (can be a scalar or NumPy array).
    b : array_like or float
        Second coordinate (can be a scalar or NumPy array).
    i : int
        Polynomial degree in the a-direction.
    j : int
        Polynomial degree in the b-direction.

    Returns
    -------
    P : ndarray or float
        The value(s) of the simplex polynomial of order (i, j) at (a, b).
    """
    # Compute the univariate Jacobi polynomials:
    #   h1 = P_i^(0,0)(a)
    #   h2 = P_j^(2*i+1,0)(b)
    h1 = JacobiP(a, 0, 0, i)
    h2 = JacobiP(b, 2*i + 1, 0, j)

    # Form the 2D simplex polynomial
    P = np.sqrt(2.0) * h1 * h2 * (1 - b) ** i
    return P

def rstoab(r, s):
    """
    Converte as coordenadas do triângulo (r, s) para (a, b).

    Parameters
    ----------
    r : array_like
        Array de coordenadas r.
    s : array_like
        Array de coordenadas s.

    Returns
    -------
    a : ndarray
        Array de coordenadas a.
    b : ndarray
        Array de coordenadas b.
    """
    r = np.asarray(r)
    s = np.asarray(s)
    a = np.empty_like(r, dtype=float)

    # For s != 1, apply the mapping; otherwise set a = -1
    mask = (s != 1)
    a[mask] = 2 * (1 + r[mask]) / (1 - s[mask]) - 1
    a[~mask] = -1

    b = s.copy()
    return a, b

def Vandermonde2D(N, r, s):
    """
    Inicializa a matriz de Vandermonde 2D,
    V2D[i,j] = phi_j(r_i, s_i) para pontos (r,s).

    Parâmetros:
    -----------
    N : int
        Grau máximo dos polinômios na base de simplex 2D.
    r : array_like, shape (M,)
        Coordenadas r dos pontos de amostragem.
    s : array_like, shape (M,)
        Coordenadas s dos pontos de amostragem.

    Retorna:
    --------
    V2D : numpy.ndarray, shape (M, (N+1)*(N+2)//2)
        Matriz de Vandermonde 2D, em que cada coluna corresponde a uma
        função base simples de graus (i,j) com i+j ≤ N.
    """

    # Inicializa a matriz V2D com zeros
    V2D = np.zeros((len(r), (N + 1) * (N + 2) // 2))

    # Transfere as coordenadas (r,s) para coordenadas (a,b) no triângulo
    a, b = rstoab(r, s)

    # Contador de coluna (índice Python começa em 0)
    sk = 0

    # Preenche cada coluna de V2D com Simplex2DP(a,b,i,j)
    # O loop percorre todas as combinações de i e j tais que i+j ≤ N
    for i in range(0, N + 1):
        for j in range(0, N + 1 - i):
            # Simplex2DP deve retornar um vetor de tamanho M com valores de φ_{i,j}(a,b)
            V2D[:, sk] = Simplex2DP(a, b, i, j)
            sk += 1

    return V2D

def Vandermonde3D(N, r, s, t):
    """
    Initialize the 3D Vandermonde matrix V3D of size (len(r), ncols),
    where ncols = (N+1)*(N+2)*(N+3)//6.

    Parameters
    ----------
    N : int
        Polynomial order.
    r, s, t : np.ndarray, shape (M,)
        Coordinates in the reference tetrahedron.

    Returns
    -------
    V3D : np.ndarray, shape (M, ncols)
        3D Vandermonde matrix, with entries V3D[i, k] = phi_k(r[i], s[i], t[i]).
    """
    # number of nodes
    M = r.size
    # number of basis functions
    ncols = (N + 1) * (N + 2) * (N + 3) // 6
    V3D = np.zeros((M, ncols), dtype=float)

    # map (r,s,t) to (a,b,c) coordinates
    a, b, c = rsttoabc(r, s, t)

    # fill Vandermonde matrix
    sk = 0
    for i in range(N + 1):
        for j in range(N + 1 - i):
            for k in range(N + 1 - i - j):
                # evaluate orthonormal polynomial on simplex
                V3D[:, sk] = Simplex3DP(a, b, c, i, j, k)
                sk += 1

    return V3D

def GradSimplex3DP(a, b, c, id, jd, kd):
    """
    Return the derivatives of the modal basis (id, jd, kd)
    on the 3D simplex at (a, b, c).

    Parameters
    ----------
    a, b, c : array_like or float
        Coordinates in the reference simplex.
    id, jd, kd : int
        Polynomial orders in the modal basis.

    Returns
    -------
    V3Dr, V3Ds, V3Dt : array_like or float
        Derivatives of the 3D Vandermonde matrix in r, s, and t.
    """
    # Evaluate Jacobi polynomials and their gradients
    fa = JacobiP(a, 0, 0, id)
    dfa = GradJacobiP(a, 0, 0, id)
    gb = JacobiP(b, 2*id + 1, 0, jd)
    dgb = GradJacobiP(b, 2*id + 1, 0, jd)
    hc = JacobiP(c, 2*(id + jd) + 2, 0, kd)
    dhc = GradJacobiP(c, 2*(id + jd) + 2, 0, kd)

    # r-derivative
    V3Dr = dfa * (gb * hc)
    if id > 0:
        V3Dr *= (0.5 * (1 - b))**(id - 1)
    if id + jd > 0:
        V3Dr *= (0.5 * (1 - c))**(id + jd - 1)

    # s-derivative
    V3Ds = 0.5 * (1 + a) * V3Dr
    tmp = dgb * (0.5 * (1 - b))**id
    if id > 0:
        tmp += -0.5 * id * (gb * (0.5 * (1 - b))**(id - 1))
    if id + jd > 0:
        tmp *= (0.5 * (1 - c))**(id + jd - 1)
    tmp = fa * (tmp * hc)
    V3Ds += tmp

    # t-derivative
    V3Dt = 0.5 * (1 + a) * V3Dr + 0.5 * (1 + b) * tmp
    tmp2 = dhc * (0.5 * (1 - c))**(id + jd)
    if id + jd > 0:
        tmp2 -= 0.5 * (id + jd) * (hc * (0.5 * (1 - c))**(id + jd - 1))
    tmp2 = fa * (gb * tmp2)
    tmp2 *= (0.5 * (1 - b))**id
    V3Dt += tmp2

    # normalize
    factor = 2**(2*id + jd + 1.5)
    V3Dr *= factor
    V3Ds *= factor
    V3Dt *= factor

    return V3Dr, V3Ds, V3Dt

def GradVandermonde3D(N, r, s, t):
    """
    Initialize the gradients of the modal basis (i,j,k) at (r,s,t) up to order N.

    Parameters
    ----------
    N : int
        Maximum polynomial order.
    r, s, t : array_like, shape (M,)
        Coordinates in the reference tetrahedron.

    Returns
    -------
    V3Dr, V3Ds, V3Dt : ndarray, shape (M, K)
        Gradient matrices in r, s, and t directions, where
        K = (N+1)*(N+2)*(N+3)/6 is the number of basis functions.
    """
    M = r.size

    # Allocate output arrays
    V3Dr = np.zeros((M, (N+1)*(N+2)*(N+3) // 6))
    V3Ds = np.zeros((M, (N+1)*(N+2)*(N+3) // 6))
    V3Dt = np.zeros((M, (N+1)*(N+2)*(N+3) // 6))

    # Map (r,s,t) to tensor-product (a,b,c) coordinates
    a, b, c = rsttoabc(r, s, t)

    # Fill in each modal gradient
    sk = 0
    for i in range(N+1):
        for j in range(N+1-i):
            for k in range(N+1-i-j):
                dr, ds, dt = GradSimplex3DP(a, b, c, i, j, k)
                V3Dr[:, sk] = dr
                V3Ds[:, sk] = ds
                V3Dt[:, sk] = dt
                sk += 1

    return V3Dr, V3Ds, V3Dt

def DMatrices3D(N, r, s, t, V):
    """
    Inicializa as matrizes de diferenciação (r,s,t) no simplex,
    avaliadas nos pontos (r,s,t) de ordem N.

    Parâmetros
    ----------
    N : int
        Ordem polinomial.
    r, s, t : array-like, shape (Np,)
        Coordenadas dos pontos de interpolação.
    V : ndarray, shape (Np, Np)
        Matriz de Vandermonde 3D de ordem N em (r,s,t).

    Retorna
    -------
    Dr, Ds, Dt : ndarray, shape (Np, Np)
        Matrizes de diferenciação em relação a r, s e t.
    """
    # Gradientes da Vandermonde no simplex
    Vr, Vs, Vt = GradVandermonde3D(N, r, s, t)
    
    # Inverte V apenas uma vez
    invV = np.linalg.inv(V)
    
    # Cálculo das matrizes de diferenciação
    Dr = Vr.dot(invV)
    Ds = Vs.dot(invV)
    Dt = Vt.dot(invV)
    
    return Dr, Ds, Dt

def Lift3D(N, R, S, T):
    """
    Compute the 3D surface-to-volume lift operator used in the DG formulation.

    Parameters
    ----------
    N : int
        Polynomial order.
    R, S, T : array_like
        Coordinates of volume nodes in reference element.

    Returns
    -------
    LIFT : ndarray, shape (Np, Nfaces*Nfp)
        The lift operator mapping surface integrals to volume.
    """
    # Initialize the face-to-volume lifting matrix
    Emat = np.zeros((Np, Nfaces * Nfp))

    # Loop over each face
    for face in range(Nfaces):
        # Extract local face coordinates r,s based on face index
        if face == 0:
            faceR = R[Fmask[:, 0]]
            faceS = S[Fmask[:, 0]]
        elif face == 1:
            faceR = R[Fmask[:, 1]]
            faceS = T[Fmask[:, 1]]
        elif face == 2:
            faceR = S[Fmask[:, 2]]
            faceS = T[Fmask[:, 2]]
        elif face == 3:
            faceR = S[Fmask[:, 3]]
            faceS = T[Fmask[:, 3]]
        else:
            raise ValueError(f"Unexpected face index: {face}")

        # Build 2D Vandermonde matrix on this face
        VFace = Vandermonde2D(N, faceR, faceS)
        # Compute face mass matrix: inv(VFace * VFace^T)
        massFace = np.linalg.inv(VFace @ VFace.T)

        # Indices in the global matrix for this face
        idr = Fmask[:, face]
        idc = np.arange(face * Nfp, (face + 1) * Nfp)

        # Accumulate into Emat
        Emat[np.ix_(idr, idc)] += massFace

    # Assemble final lift operator: V * (V^T * Emat)
    LIFT = V @ (V.T @ Emat)
    return LIFT

def GradH3D(U):
    """
    Compute local elemental physical spatial derivatives of U in 3D.

    Parameters
    ----------
    U : ndarray, shape (Np,)
        Valores de U nos nós do elemento.

    Returns
    -------
    dUdx, dUdy, dUdz : ndarray, shape (Np,)
        Derivadas físicas espaciais de U em x, y e z.
    """
    # derivadas no tetraedro de referência
    dUdr = Dr @ U
    dUds = Ds @ U
    dUdt = Dt @ U

    # aplicação da regra da cadeia para derivadas físicas
    dUdx = rx * dUdr + sx * dUds + tx * dUdt
    dUdy = ry * dUdr + sy * dUds + ty * dUdt
    dUdz = rz * dUdr + sz * dUds + tz * dUdt

    return dUdx, dUdy, dUdz

def DivH3D(Ux, Uy, Uz):
    """
    Compute local elemental physical spatial divergence of (Ux, Uy, Uz) in 3D.

    Parameters
    ----------
    Ux, Uy, Uz : ndarray, shape (Np,)
        Components of the vector field U at the element's nodes.

    Returns
    -------
    divU : ndarray, shape (Np,)
        The divergence ∇·U at each node of the element.
    """
    # dUx/dx
    ddr = Dr @ Ux
    dds = Ds @ Ux
    ddt = Dt @ Ux
    divU = rx * ddr + sx * dds + tx * ddt

    # dUy/dy
    ddr = Dr @ Uy
    dds = Ds @ Uy
    ddt = Dt @ Uy
    divU += ry * ddr + sy * dds + ty * ddt

    # dUz/dz
    ddr = Dr @ Uz
    dds = Ds @ Uz
    ddt = Dt @ Uz
    divU += rz * ddr + sz * dds + tz * ddt

    return divU

def CurlH3D(Ux, Uy, Uz):
    """
    Compute local elemental physical spatial curl of (Ux, Uy, Uz) in 3D.

    Parameters
    ----------
    Ux, Uy, Uz : ndarray, shape (Np,)
        Components of the vector field U at the element's nodes.

    Returns
    -------
    curlx, curly, curlz : ndarray, shape (Np,)
        The components of curl U at each node of the element.
    """
    # ∂U_x on reference tetrahedron
    ddr = Dr @ Ux
    dds = Ds @ Ux
    ddt = Dt @ Ux
    curly =  rz * ddr + sz * dds + tz * ddt
    curlz = - (ry * ddr + sy * dds + ty * ddt)

    # ∂U_y on reference tetrahedron
    ddr = Dr @ Uy
    dds = Ds @ Uy
    ddt = Dt @ Uy
    curlx = - (rz * ddr + sz * dds + tz * ddt)
    curlz +=   rx * ddr + sx * dds + tx * ddt

    # ∂U_z on reference tetrahedron
    ddr = Dr @ Uz
    dds = Ds @ Uz
    ddt = Dt @ Uz
    curlx +=   ry * ddr + sy * dds + ty * ddt
    curly -=   rx * ddr + sx * dds + tx * ddt

    return curlx, curly, curlz

def generate_mesh(refinement_level, a, b, c, L, Ls):

    #Dimensões do domínio físico (sem PML)
    xf = a - L - Ls
    yf = b - L - Ls
    zf = c - L - Ls
    #Dimensões da região TFSF
    xS = a - L
    yS = b - L
    zS = c - L

    #Definindo os vértices do cubo
    points =[
        (-a, -b, -c), #Ponto 0
        (a, -b, -c),  #Ponto 1
        (a, b, -c),   #Ponto 2
        (-a, b, -c),  #Ponto 3
        (-a, -b, c),  #Ponto 4
        (a, -b, c),   #Ponto 5
        (a, b, c),    #Ponto 6
        (-a, b, c),   #Ponto 7, até aqui todos os pontos são externos (fora da PML)
        (-xf, -yf, -zf), #Ponto 8
        (xf, -yf, -zf),  #Ponto 9
        (xf, yf, -zf),   #Ponto 10
        (-xf, yf, -zf),  #Ponto 11
        (-xf, -yf, zf),  #Ponto 12
        (xf, -yf, zf),   #Ponto 13
        (xf, yf, zf),    #Ponto 14
        (-xf, yf, zf),   #Ponto 15
        (-xS, -yS, -zS), #Ponto 16
        (xS, -yS, -zS),  #Ponto 17
        (xS, yS, -zS),   #Ponto 18
        (-xS, yS, -zS),  #Ponto 19
        (-xS, -yS, zS),  #Ponto 20
        (xS, -yS, zS),   #Ponto 21
        (xS, yS, zS),    #Ponto 22
        (-xS, yS, zS)    #Ponto 23
    ]

    facets = [
        # Base inferior
        [[0, 1, 2], [0, 2, 3]],
        # Base superior
        [[4, 5, 6], [4, 6, 7]],
        # Faces laterais
        [[0, 1, 5], [0, 5, 4]],
        [[1, 2, 6], [1, 6, 5]],
        [[2, 3, 7], [2, 7, 6]],
        [[3, 0, 4], [3, 4, 7]],
        #Cubo médio
        [[16, 17, 18], [16, 18, 19]],
        [[20, 21, 22], [20, 22, 23]],
        [[16, 17, 21], [16, 21, 20]],
        [[17, 18, 22], [17, 22, 21]],
        [[18, 19, 23], [18, 23, 22]],
        [[19, 16, 20], [19, 20, 23]],
        #Cubo menor 
        [[8, 9, 10], [8, 10, 11]],
        [[12, 13, 14], [12, 14, 15]],
        [[8, 9, 13], [8, 13, 12]],
        [[9, 10, 14], [9, 14, 13]],
        [[10, 11, 15], [10, 15, 14]],
        [[11, 8, 12], [11, 12, 15]]
    ]

    mesh_info = MeshInfo()
    mesh_info.set_points(points)
    # As facetas devem ser uma lista plana de listas ou tuplas
    flattened_facets = [facet for pair in facets for facet in pair]
    mesh_info.set_facets(flattened_facets)

    mesh = build(mesh_info, max_volume = 1 / (refinement_level**2))
    
    vertices = np.array(mesh.points)
    elements = np.array(mesh.elements)

    """print(f"Nós: {np.array(vertices)}")
    print(f'Elementos: {np.array(elements)}')"""

    sorted_elements = np.sort(elements)

    edge_dict = {}  # Dicionário para armazenar arestas únicas e seus índices
    edge_list = []  # Lista para armazenar as arestas globais na ordem de criação
    edge_map = []  # Mapeamento local de arestas para elementos

    current_edge_index = 0

    for element in sorted_elements:
        # Define as arestas do triângulo com base na convenção específica
        local_edges = [
            (element[0], element[1]),  # Aresta 0: nó 0 -> nó 1
            (element[0], element[2]),  # Aresta 1: nó 0 -> nó 2
            (element[0], element[3]),  # Aresta 2: nó 0-> nó 3
            (element[1], element[2]),  # Aresta 3: nó 1 -> nó 2
            (element[1], element[3]),  # Aresta 4: nó 1 -> nó 3
            (element[2], element[3])   # Aresta 5: nó 2 -> nó 3
        ]

        # Normaliza as arestas (min, max) para garantir consistência
        local_edges = [tuple(edge) for edge in local_edges]

        # Mapeia as arestas para os índices globais
        local_edge_indices = []
        for edge in local_edges:
            if edge not in edge_dict:
                edge_dict[edge] = current_edge_index
                edge_list.append(edge)
                current_edge_index += 1
            local_edge_indices.append(edge_dict[edge])

        edge_map.append(local_edge_indices)

    # Converte edge_list e edge_map para arrays numpy
    edges = np.array(edge_list)
    edge_map = np.array(edge_map)


    return vertices, elements, edges, edge_map

def plot_thetaedral_mesh(vertices, elements, edges, edge_map):
    fig = plt.figure(figsize=(10, 10))
    ax = fig.add_subplot(111, projection='3d')

    #Plotar os elementos tetraédricos
    for elem_idx, elem in enumerate(elements):
        elem_vertices = vertices[elem, :] #Pegando os vértices do tetraedro
        faces = [elem_vertices[[0, 1, 2]],
                 elem_vertices[[0, 1, 3]],
                 elem_vertices[[0, 2, 3]],
                 elem_vertices[[1, 2, 3]]]
        
        poly3d = Poly3DCollection(faces, edgecolor='black', alpha=0.3)
        ax.add_collection3d(poly3d)

        #Adicionar as arestas dos elementos
        local_edges = edge_map[elem_idx]
        for edge_local_idx, edge_global_idx in enumerate(local_edges):
            edge = edges[edge_global_idx]
            start_vertex, end_vertex = vertices[edge[0]], vertices[edge[1]]

            ax.plot(
                [start_vertex[0], end_vertex[0]],
                [start_vertex[1], end_vertex[1]],
                [start_vertex[2], end_vertex[2]],
                color='blue', linewidth=0.8
            )

            #Adiciona rótulo na aresta
            mid_x = (start_vertex[0] + end_vertex[0]) / 2
            mid_y = (start_vertex[1] + end_vertex[1]) / 2
            mid_z = (start_vertex[2] + end_vertex[2]) / 2
            ax.text(mid_x, mid_y, mid_z, f'{edge_global_idx}', fontsize=8, color='red', ha='center')
    
    #Plotar os nós
    ax.scatter(vertices[:, 0], vertices[:, 1], vertices[:, 2], c='red', s=20, label = 'Nós')

    #Anotar os índices dos nós
    for i, (x, y, z) in enumerate(vertices):
        ax.text(x, y, z, f'{i}', fontsize=8, ha='right', va='top')

    ax.set_title('Malha de Elementos Finitos')
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('z')
    plt.legend

def export_mesh(vertices, elements):

    nodes = vertices

    # Criando a estrutura da malha no VTK
    points = vtk.vtkPoints()
    for node in nodes:
        points.InsertNextPoint(node)

    # Criando células da malha (tetraedros)
    unstructuredGrid = vtk.vtkUnstructuredGrid()
    unstructuredGrid.SetPoints(points)

    for elem in elements:
        tetra = vtk.vtkTetra()
        for i in range(4):
            tetra.GetPointIds().SetId(i, elem[i])
        unstructuredGrid.InsertNextCell(tetra.GetCellType(), tetra.GetPointIds())

    # Criando o writer para salvar no formato VTU
    writer = vtk.vtkXMLUnstructuredGridWriter()
    writer.SetFileName("malhanova.vtu")
    writer.SetInputData(unstructuredGrid)
    writer.Write()

    print("Arquivo VTU salvo como 'malha.vtu'.")

def get_EToV(dim, index_based=1):
    """
    # ### ``gmsh/model/mesh/getElements``
    # 
    # > Get the elements classified on the entity of dimension `dim` and tag ``tag``. If ``tag`` < 0, 
    #   get the elements for all entities of dimension ``dim``. If ``dim`` and ``tag`` are negative, 
    #   get all the elements in the mesh. ``elementTypes`` contains the MSH types of the elements 
    #   (e.g. ``2`` for 3-node triangles: see ``getElementProperties`` to obtain the properties for 
    #   a given element type). 
    # 
    #   ``elementTags`` is a vector of the same length as ``elementTypes``; 
    #   each entry is a vector containing the tags (unique, strictly positive identifiers) of the 
    #   elements of the corresponding type. ``nodeTags`` is also a vector of the same length as 
    #   ``elementTypes``; each entry is a vector of length equal to the number of elements of the 
    #   given type times the number N of nodes for this type of element, that contains the node tags 
    #   of all the elements of the given type, concatenated: [e1n1, e1n2, ..., e1nN, e2n1, ...].
    # 
    # > Input: ``dim`` = -1 (integer), ``tag`` = -1 (integer)
    # > Output: ``elementTypes`` (vector of integers), ``elementTags`` (vector of vectors of sizes), 
    #           ``nodeTags`` (vector of vectors of sizes)
    # > Return: -
    # 
    # ### `gmsh/model/mesh/getElementProperties`
    # 
    # > Get the properties of an element of type `elementType`: its name (``elementName``), 
    #   dimension (``dim``), order (``order``), number of nodes (``numNodes``), local coordinates 
    #   of the nodes in the reference element (``localNodeCoord`` vector, of length ``dim`` times 
    #   ``numNodes``) and number of primary (first order) nodes (``numPrimaryNodes``).
    # 
    # > Input: ``elementType`` (integer)  
    # > Output: ``elementName`` (string), ``dim`` (integer), ``order`` (integer), ``numNodes`` (integer), 
    #           ``localNodeCoord`` (vector of doubles), ``numPrimaryNodes`` (integer)  
    # > Return: -
    """
    # Obter os elementos da malha
    elemTypes, elemTags, elemNodeTags = gmsh.model.mesh.getElements(dim)
    EToV = []

    for elemType, elemNode in zip(elemTypes, elemNodeTags):
        # Obter as propriedades do elemento
        _, _, _, nodes_per_element, _, _ = gmsh.model.mesh.getElementProperties(elemType)
        numElements = len(elemNode) // nodes_per_element
        
        # Para cada elemento, extrair os índices dos vértices e converter para int
        for i in range(numElements):
            conn_i = elemNode[nodes_per_element * i: nodes_per_element * (i + 1)]
            EToV.append([int(tag) for tag in conn_i])

    # Se index_based for False, converter os índices para base 0
    if not index_based:
        EToV = [np.array(conn) - 1 for conn in EToV]

    # Retornar a matriz de conectividade. Dim: K x N
    # Kx4 para tetraedros; Kx3 para triângulos)
    return np.array(EToV, dtype=int)  

def extract_VX_VY_VZ(nodes_data):
    # Ordena os nós pelo tag (ID do nó)
    sorted_nodes = sorted(nodes_data.items())

    # Extrai coordenadas ordenadas
    VX = np.array([data["xg"][0] for _, data in sorted_nodes])
    VY = np.array([data["xg"][1] for _, data in sorted_nodes])
    VZ = np.array([data["xg"][2] for _, data in sorted_nodes])

    return VX, VY, VZ

def get_nodes_data(dim, BOUNDARY=None, INTERFACES=None):
    # 1. Dicionário para Mapeamento inicial de nós:
    NodeTags, NodeCoords, _ = gmsh.model.mesh.getNodes()

    # 2. Converte os NodeTags em inteiros nativos
    NodeTags = NodeTags.astype(np.uint64)
    NodeTags = [int(tag) for tag in NodeTags]

    # 2.1 Todos os nós começam com a condição de contorno "Free" e valor None.
    node_bc_map = {
        tag:{
            'tag': None,
            'type': 'Free',
            'value': None,
            'name': 'free_node'
        }
        for tag in NodeTags
    }

    # 3. Atualização do mapeamento de nós:
    # 3.1 Para cada condição de contorno em BOUNDARY, os nós associados ao grupo físico correspondente 
    # são atualizados com o tipo (type) e valor (value) dessa condição.
    if BOUNDARY is not None:
        for bc in BOUNDARY:
            # Obtenha os nós associados ao grupo físico especificado no bc
            condition_NodeTags, _ = gmsh.model.mesh.getNodesForPhysicalGroup(dim, tag=bc['tag'])
            
            # Atualiza o mapeamento de nós com a condição de contorno correspondente
            for node in condition_NodeTags:
                node_bc_map[node] = {
                    'tag': bc['tag'],
                    'type': bc['type'],
                    'value': bc['value'],
                    'name': bc['name']
                }
    
    # 3.2 Para cada condição de contorno em INTERFACES, os nós associados ao grupo físico correspondente 
    # são atualizados com o tipo (type) e valor (value) dessa condição.
    if INTERFACES is not None:
        for bc in INTERFACES:
            # Obtenha os nós associados ao grupo físico especificado no bc
            condition_NodeTags, _ = gmsh.model.mesh.getNodesForPhysicalGroup(dim, tag=bc['tag'])
            
            # Atualiza o mapeamento de nós com a condição de contorno correspondente
            for node in condition_NodeTags:
                node_bc_map[node] = {
                    'tag': bc['tag'],
                    'type': bc['type'],
                    'value': bc['value'],
                    'name': bc['name']
                }

    # 4. Estrutura final:
    # A lista dict_nodes contém informações completas sobre cada nó, incluindo suas coordenadas 
    # globais e as condições de contorno associadas.
    dict_nodes = {
        int(node): {
            "xg": (
                float(NodeCoords[3*i]),
                float(NodeCoords[3*i + 1]),
                float(NodeCoords[3*i + 2])
            ),
            "bc": node_bc_map[int(node)]
        } 
        for i, node in enumerate(NodeTags)
    }

    return dict_nodes

def get_faces(tet):
    return [
        [tet[0], tet[1], tet[2]],
        [tet[0], tet[1], tet[3]],
        [tet[0], tet[2], tet[3]],
        [tet[1], tet[2], tet[3]]
    ]

def get_edges(EToV):
    """
    Retorna arestas externas únicas da malha convexa formada por EToV.

    Parâmetros:
        EToV : ndarray (K x 4) — conectividade dos tetraedros
        vx, vy, vz : ndarray — coordenadas dos vértices
        atol : float — tolerância para evitar duplicação numérica

    Retorna:
        edges : list of tuples (i, j) — arestas únicas da geometria externa
    """
    from collections import defaultdict
    edge_dict = defaultdict(int)
    edge_set = set()

    # Cada tetraedro tem 6 arestas: todas combinações de 2 entre 4 vértices
    for tet in EToV:
        for i in range(4):
            for j in range(i+1, 4):
                ni, nj = sorted((tet[i], tet[j]))
                edge_dict[(ni, nj)] += 1
                edge_set.add((ni, nj))

    # Arestas que aparecem só uma vez são externas (na superfície da malha convexa)
    # outer_edges = [edge for edge, count in edge_dict.items() if count == 1]

    return sorted(edge_set)

def plot_tetrahedral_mesh2(VX, VY, VZ, EToV):
    # Arestas do cubo (pares de índices dos nós)
    cube_edges = get_edges(EToV)

    fig = plt.figure(figsize=(10, 8))
    colors = plt.cm.tab10(np.linspace(0, 1, len(EToV)))
    for k, elem in enumerate(EToV):
        ax = fig.add_subplot(2, 3, k + 1, projection='3d')
        # ax.view_init(elev=0, azim=30)

        # 1. Desenhar arestas do cubo (tracejadas)
        for i, j in cube_edges:
            x = [VX[i], VX[j]]
            y = [VY[i], VY[j]]
            z = [VZ[i], VZ[j]]
            ax.plot(x, y, z, 'k--', linewidth=0.8, zorder=0)

        # 2. Desenhar tetraedro colorido
        faces = get_faces(elem)
        for face in faces:
            verts = [[VX[i], VY[i], VZ[i]] for i in face]
            poly = Poly3DCollection([verts], facecolor=colors[k], edgecolor='k', alpha=0.7)
            ax.add_collection3d(poly)

        # 3. Adicionar vértices e rótulos
        ax.scatter(VX[elem], VY[elem], VZ[elem], color='red')
        for idx in elem:
            ax.text(VX[idx]+0.05, VY[idx]+0.05, VZ[idx]+0.05, str(idx), color='black', fontsize=10)

        # 4. Rótulo do ID do elemento
        centroid = np.mean([[VX[i], VY[i], VZ[i]] for i in elem], axis=0)
        ax.text(*centroid, f't{k}', color='black', ha='center', fontsize=11)

        # 5. Ajuste dos eixos
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
        ax.set_zlim(0, 1)
        ax.set_box_aspect([1, 1, 1])
        ax.axis('off')

    fig.suptitle("Unit cube", fontsize=12)
    plt.tight_layout()

def GeometricFactors3D(x, y, z, Dr, Ds, Dt):
    """
    Compute the metric elements for the local mappings of 3D elements.

    Parameters
    ----------
    x, y, z : ndarray, shape (Np,)
        Coordinates of nodes in physical space.
    Dr, Ds, Dt : ndarray, shape (Np, Np)
        Differentiation matrices in the reference element directions r, s, t.

    Returns
    -------
    rx, sx, tx, ry, sy, ty, rz, sz, tz, J : ndarrays, shape (Np,)
        Geometric factors and Jacobian at each node.
    """
    # Compute physical derivatives
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

def Normals3D():
    """
    Compute outward pointing normals at element faces and surface Jacobians.

    Returns
    -------
    nx, ny, nz : ndarray of shape (4*Nfp, K)
        Unit outward normal components on each face node.
    sJ : ndarray of shape (4*Nfp, K)
        Surface Jacobian (scaled by volume Jacobian at face nodes).
    """
    # Compute geometric factors
    rx, sx, tx, ry, sy, ty, rz, sz, tz, J = GeometricFactors3D(x, y, z, Dr, Ds, Dt)

    # Interpolate geometric factors to face nodes
    face_idx = Fmask.flatten('F')
    frx = rx[face_idx, :]
    fsx = sx[face_idx, :]
    ftx = tx[face_idx, :]
    fry = ry[face_idx, :]
    fsy = sy[face_idx, :]
    fty = ty[face_idx, :]
    frz = rz[face_idx, :]
    fsz = sz[face_idx, :]
    ftz = tz[face_idx, :]

    # Initialize normal and surface Jacobian arrays
    totalF = 4 * Nfp
    nx = np.zeros((totalF, K))
    ny = np.zeros((totalF, K))
    nz = np.zeros((totalF, K))

    # Face index ranges
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

    # Normalize normals and compute surface Jacobian
    sJ = np.sqrt(nx**2 + ny**2 + nz**2)
    nx /= sJ
    ny /= sJ
    nz /= sJ

    # Scale sJ by volume Jacobian at face nodes
    sJ = sJ * J[face_idx, :]

    return nx, ny, nz, sJ

def tiConnect3D(EToV):
    """
    Tetrahedral face connectivity algorithm (Toby Isaac).

    Parameters
    ----------
    EToV : array_like, shape (K, 4)
        Element-to-vertex mapping. Node indices devem ser 0-based.

    Returns
    -------
    EToE : ndarray, shape (K, 4)
        Conectividade elemento-elemento. EToE[e, f] é o índice do elemento
        vizinho ao elemento e na face f (0-based).
    EToF : ndarray, shape (K, 4)
        Conectividade elemento-face. EToF[e, f] é o número da face local no
        elemento vizinho que faz contato na face f do elemento e (0-based).
    """
    Nnodes = VX.shape[0]
                
    # create list of all faces 1, then 2, 3 & 4
    fnodes = np.concatenate(
            (EToV[:,[0, 1, 2]],    # face oposta ao vértice 3
             EToV[:,[0, 1, 3]],    # face oposta ao vértice 2
             EToV[:,[1, 2, 3]],    # face oposta ao vértice 0  
             EToV[:,[0, 2, 3]])    # face oposta ao vértice 1
        )

    # ordena vértices de cada face para identificação única
    fnodes = np.sort(fnodes, axis=1)

    # set up default element to element and Element to faces connectivity
    EToE = np.outer(np.arange(K, dtype=int), np.ones((Nfaces,), dtype=int))
    EToF = np.outer(np.ones((K,), dtype=int), np.arange(Nfaces, dtype=int))

    # Gera identificadores únicos por face
    ids = fnodes[:, 0] * Nnodes * Nnodes + fnodes[:, 1] * Nnodes + fnodes[:, 2]

    spNodeToNode = np.concatenate([
        ids.reshape(-1, 1),
        np.arange(Nfaces * K).reshape(-1, 1),
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

def BuildMaps3D():

    """
    Esta função constrói os mapeamentos entre os nós das faces de elementos. 

    Retorna:
        mapM   : array inteiro de tamanho (K*Nfp*Nfaces,), mapeamento “mais interno” das faces
        mapP   : array inteiro de tamanho (K*Nfp*Nfaces,), mapeamento “par” correspondente
        vmapM  : array inteiro de tamanho (K*Nfp*Nfaces,), posição de cada nó de face “mais interno” no volume
        vmapP  : array inteiro de tamanho (K*Nfp*Nfaces,), posição de cada nó de face “par” no volume
        vmapB  : array inteiro contendo os índices de vmapM que são fronteira (vmapP == vmapM)
        mapB   : array inteiro contendo os índices lineares das faces de fronteira (mapM)
    """

    # Cria os índices globais dos nós em cada elemento (nodeids tem shape (Np, K))
    nodeids = np.linspace(0, K * Np - 1, K * Np, dtype=int).reshape((K, Np)).T

    # Inicializa os arrays que vão guardar os índices globais dos nós nas faces
    vmapM = np.zeros((Nfp, Nfaces, K), dtype=int)
    vmapP = np.zeros((Nfp, Nfaces, K), dtype=int)

    # Cria o mapeamento base mapM e inicializa mapP com o mesmo formato
    n = K * Nfp * Nfaces
    mapM = np.linspace(0, n-1, n, dtype=int)[:, None]
    mapP = np.copy(mapM).reshape((Nfp, Nfaces, K))  # será atualizado com os vizinhos

    # Preenche vmapM com os índices globais dos nós nas faces de cada elemento
    for k1 in range(K):
        for f1 in range(Nfaces):
            vmapM[:, f1, k1] = nodeids[Fmask[:, f1], k1]

    # Matriz auxiliar de uns, usada para broadcasting nas distâncias
    one = np.ones((1, Nfp))

    # Loop sobre cada elemento e cada face, identificando vizinhos e mapeando os nós entre faces opostas
    for k1 in range(K):
        for f1 in range(Nfaces):
            k2 = EToE[k1, f1]  # elemento vizinho
            f2 = EToF[k1, f1]  # face correspondente no vizinho

            # Índices globais dos nós da face atual (menor lado)
            vidM = vmapM[:, f1, k1]
            # Índices globais dos nós da face do vizinho (maior lado)
            vidP = vmapM[:, f2, k2]

            # Coordenadas dos nós das duas faces
            x1 = x.flatten('F')[vidM]
            y1 = y.flatten('F')[vidM]
            z1 = z.flatten('F')[vidM]
            x2 = x.flatten('F')[vidP]
            y2 = y.flatten('F')[vidP]
            z2 = z.flatten('F')[vidP]
            x1, y1, z1 = x1[:, None] @ one, y1[:, None] @ one, z1[:, None] @ one
            x2, y2, z2 = x2[:, None] @ one, y2[:, None] @ one, z2[:, None] @ one

            # Matriz de distâncias quadradas entre nós das duas faces
            D = (x1 - x2.T)**2 + (y1 - y2.T)**2 + (z1 - z2.T)**2

            # Identifica pares de nós que coincidem (dentro da tolerância)
            idM, idP = np.where(np.sqrt(np.abs(D)) < NODETOL)

            # Atualiza os mapeamentos com os índices corretos dos vizinhos
            vmapP[idM, f1, k1] = vidP[idP]
            mapP[idM, f1, k1] = idP + f2 * Nfp + k2 * Nfaces * Nfp

    # Achata os arrays tridimensionais para vetores coluna
    vmapM_flat = vmapM.flatten(order='F')
    vmapP_flat = vmapP.flatten(order='F')
    mapM = mapM.flatten(order='F')
    mapP_flat = mapP.flatten(order='F')

    # Identifica os nós de fronteira: onde o vizinho é o próprio nó (i.e., vmapP == vmapM)
    mapB = np.where(vmapP_flat == vmapM_flat)[0]
    vmapB = vmapM_flat[mapB]

    return mapM, mapP_flat, vmapM_flat, vmapP_flat, vmapB, mapB

def rungekutta():
    """
    Calcula os coeficientes do runge-kutta de 4ª ordem
    """

    rk4a = np.array([
    0.0,
    -567301805773.0 / 1357537059087.0,
    -2404267990393.0 / 2016746695238.0,
    -3550918686646.0 / 2091501179385.0,
    -1275806237668.0 / 842570457699.0
    ])

    rk4b = np.array([
        1432997174477.0 / 9575080441755.0,
        5161836677717.0 / 13612068292357.0,
        1720146321549.0 / 2090206949498.0,
        3134564353537.0 / 4481467310338.0,
        2277821191437.0 / 14882151754819.0
    ])

    rk4c = np.array([
        0.0,
        1432997174477.0 / 9575080441755.0,
        2526269341429.0 / 6820363962896.0,
        2006345519317.0 / 3224310063776.0,
        2802321613138.0 / 2924317926251.0
    ])

    return rk4a, rk4b, rk4c

def MaxwellRHS3D(Hx, Hy, Hz, Ex, Ey, Ez, PHx, PHy, PHz, PEx, PEy, PEz, timelocal):
    """
    Avalia o lado direito da equação de Maxwell no caso 3D.

    Parameters:
    -----------
    Hx, Hy, Hz, Ex, Ey, Ez : ndarray
        Field arrays of shape (Np, K)
    alpha : float, optional
        Upwind parameter (default is 1.0)

    Returns:
    --------
    rhsHx, rhsHy, rhsHz, rhsEx, rhsEy, rhsEz : ndarray
        Right-hand side arrays of shape (Np, K)
    """

    Hx_flat = Hx.flatten(order='F')
    Hy_flat = Hy.flatten(order='F')
    Hz_flat = Hz.flatten(order='F')
    Ex_flat = Ex.flatten(order='F')
    Ey_flat = Ey.flatten(order='F')
    Ez_flat = Ez.flatten(order='F')

    # Flatten normals and scales (assumed shape Nfp*Nfaces x K)
    nx_flat = nx.flatten(order='F')
    ny_flat = ny.flatten(order='F')
    nz_flat = nz.flatten(order='F')

    # Form field differences at faces
    dHx = Hx_flat[vmapM] - Hx_flat[vmapP]
    dHy = Hy_flat[vmapM] - Hy_flat[vmapP]
    dHz = Hz_flat[vmapM] - Hz_flat[vmapP]
    dEx = Ex_flat[vmapM] - Ex_flat[vmapP]
    dEy = Ey_flat[vmapM] - Ey_flat[vmapP]
    dEz = Ez_flat[vmapM] - Ez_flat[vmapP]

    '''# Apply reflective boundary conditions (Ez+ = -Ez-)
    dHx[mapB] = 0
    dHy[mapB] = 0
    dHz[mapB] = 0
    dEx[mapB] = 2 * Ex_flat[vmapB]
    dEy[mapB] = 2 * Ey_flat[vmapB]
    dEz[mapB] = 2 * Ez_flat[vmapB]'''

    # Upwinding parameter
    alpha = 1

    #Impondo a fonte por TFSF
    if timelocal < FinalTime:
        #Definindo a onda incidente com o sinal dependendo da região:
        Eyinc = np.cos(omega * timelocal - k * x.flatten('F'))
        Hzinc = np.cos(omega * timelocal - k * x.flatten('F'))

        '''dEy[maskScatterRegion] += Eyinc[vmapP[maskScatterRegion]]
        dEy[maskTotalregion] -= Eyinc[vmapP[maskTotalregion]]
        dHz[maskScatterRegion] += Hzinc[vmapP[maskScatterRegion]]
        dHz[maskTotalregion] -= Hzinc[vmapP[maskTotalregion]]'''
        
        dEy[maskScatterRegion] = Ey_flat[vmapM[maskScatterRegion]] - (Ey_flat[vmapP[maskScatterRegion]] - Eyinc[vmapP[maskScatterRegion]])
        dHz[maskScatterRegion] = Hz_flat[vmapM[maskScatterRegion]] - (Hz_flat[vmapP[maskScatterRegion]] - Hzinc[vmapP[maskScatterRegion]])
        
        dEy[maskTotalRegion] = Ey_flat[vmapM[maskTotalRegion]] - (Ey_flat[vmapP[maskTotalRegion]] + Eyinc[vmapP[maskTotalRegion]])
        dHz[maskTotalRegion] = Hz_flat[vmapM[maskTotalRegion]] - (Hz_flat[vmapP[maskTotalRegion]] + Hzinc[vmapP[maskTotalRegion]])


    # Compute normal dot differences
    ndotdH = nx_flat * dHx + ny_flat * dHy + nz_flat * dHz
    ndotdE = nx_flat * dEx + ny_flat * dEy + nz_flat * dEz

    # Compute numerical fluxes
    fluxHx = ny_flat * dEz - nz_flat * dEy - alpha * (dHx - ndotdH * nx_flat)
    fluxHy = nz_flat * dEx - nx_flat * dEz - alpha * (dHy - ndotdH * ny_flat)
    fluxHz = nx_flat * dEy - ny_flat * dEx - alpha * (dHz - ndotdH * nz_flat)
    fluxEx = -ny_flat * dHz + nz_flat * dHy - alpha * (dEx - ndotdE * nx_flat)
    fluxEy = -nz_flat * dHx + nx_flat * dHz - alpha * (dEy - ndotdE * ny_flat)
    fluxEz = -nx_flat * dHy + ny_flat * dHx - alpha * (dEz - ndotdE * nz_flat)

    # Reshape fields back to (Np, K) for curl operations
    # Note: Curl3D expects (Np, K) arrays
    fluxHx = fluxHx.reshape((Nfp * Nfaces, K), order = 'F')
    fluxHy = fluxHy.reshape((Nfp * Nfaces, K), order = 'F')
    fluxHz = fluxHz.reshape((Nfp * Nfaces, K), order = 'F')
    fluxEx = fluxEx.reshape((Nfp * Nfaces, K), order = 'F')
    fluxEy = fluxEy.reshape((Nfp * Nfaces, K), order = 'F')
    fluxEz = fluxEz.reshape((Nfp * Nfaces, K), order = 'F')

    # Evaluate local spatial derivatives
    curlHx, curlHy, curlHz = CurlH3D(Hx, Hy, Hz)
    curlEx, curlEy, curlEz = CurlH3D(Ex, Ey, Ez)

    # Assemble RHS with lift operator
    # Fscale and LIFT are global arrays (faces-scalings and lifting matrix) #MUDAR PRA INSERIR A PML
    rhsHx = -curlEx + LIFT @ ((Fscale * fluxHx) / 2) + PHx - (sigma_y + sigma_z - sigma_x) * Hx
    rhsHy = -curlEy + LIFT @ ((Fscale * fluxHy) / 2) + PHy - (sigma_x + sigma_z - sigma_y) * Hy
    rhsHz = -curlEz + LIFT @ ((Fscale * fluxHz) / 2) + PHz - (sigma_x + sigma_y - sigma_z) * Hz
    rhsEx =  curlHx + LIFT @ ((Fscale * fluxEx) / 2) - PEx - (sigma_y + sigma_z - sigma_x) * Ex
    rhsEy =  curlHy + LIFT @ ((Fscale * fluxEy) / 2) - PEy - (sigma_x + sigma_z - sigma_y) * Ey
    rhsEz =  curlHz + LIFT @ ((Fscale * fluxEz) / 2) - PEz - (sigma_x + sigma_y - sigma_z) * Ez
    #Para as variáveis auxiliares
    rhsPHx = -sigma_x * PHx + (sigma_x * sigma_y + sigma_x * sigma_z - sigma_y*sigma_z - sigma_x**2) * Hx
    rhsPHy = -sigma_y * PHy + (sigma_x * sigma_y + sigma_y * sigma_z - sigma_x*sigma_z - sigma_y**2) * Hy
    rhsPHz = -sigma_z * PHz + (sigma_x * sigma_z + sigma_y * sigma_z - sigma_x*sigma_y - sigma_z**2) * Hz
    rhsPEx = -sigma_x * PEx - (sigma_x * sigma_y + sigma_x * sigma_z - sigma_y*sigma_z - sigma_x**2) * Ex
    rhsPEy = -sigma_y * PEy - (sigma_x * sigma_y + sigma_y * sigma_z - sigma_x*sigma_z - sigma_y**2) * Ey
    rhsPEz = -sigma_z * PEz - (sigma_x * sigma_z + sigma_y * sigma_z - sigma_x*sigma_y - sigma_z**2) * Ez

    return rhsHx, rhsHy, rhsHz, rhsEx, rhsEy, rhsEz, rhsPHx, rhsPHy, rhsPHz, rhsPEx, rhsPEy, rhsPEz

def dtscale3D():
    """
    Faz o escalonamento temporal para fins de estabilidade
    """

    vmask1 = np.where(abs(s + r  + t + 2) < NODETOL)[0]
    vmask2 = np.where(abs(r - 1) < NODETOL)[0]
    vmask3 = np.where(abs(s -1) < NODETOL)[0]
    vmask4 = np.where(abs(t -1) < NODETOL)[0]
    vmask = np.concatenate([vmask1, vmask2, vmask3, vmask4])
    vx, vy, vz = x[vmask], y[vmask], z[vmask]

    #Calculando semi-perímetro e área
    len1 = np.sqrt((vx[0] - vx[1])**2 + (vy[0] - vy[1])**2 + (vz[0] - vz[1])**2)
    len2 = np.sqrt((vx[1] - vx[2])**2 + (vy[1] - vy[2])**2 + (vz[1] - vz[2])**2)
    len3 = np.sqrt((vx[2] - vx[0])**2 + (vy[2] - vy[0])**2 + (vz[2] - vz[0])**2)
    sper = (len1 + len2 + len3)/2
    Area = np.sqrt(sper * (sper-len1) * (sper-len2) * (sper-len3))
    #Calcula a escala usando o raio do círculo inscrito
    dtscale = Area/sper

    return dtscale

def Maxwell3D(Hx, Hy, Hz, Ex, Ey, Ez, PHx, PHy, PHz, PEx, PEy, PEz, FinalTime):
    """
    Integra o modo TM das equações de Maxwell até o tempo FinalTime
    começando com as condições iniciais Hx, Hy e Ez.

    Parâmetros:
    -----------
    Hx, Hy, Ez : np.ndarray
        Campos magnéticos nos eixos x e y e campo elétrico no eixo z.
    FinalTime : float
        Tempo final de integração.

    Retorna:
    --------
    Hx, Hy, Ez : np.ndarray
        Campos atualizados até o tempo FinalTime.
    time : float
        Tempo final (deverá ser igual a FinalTime).
    """
    # Tempo inicial
    time = 0.0

    # Inicializa armazenamento dos resíduos do método de Runge-Kutta
    resHx = np.zeros((Np, K))
    resPHx = np.zeros((Np, K))
    resHy = np.zeros((Np, K))
    resPHy = np.zeros((Np, K))
    resHz = np.zeros((Np, K))
    resPHz = np.zeros((Np, K))
    resEx = np.zeros((Np, K))
    resPEx = np.zeros((Np, K))
    resEy = np.zeros((Np, K))
    resPEy = np.zeros((Np, K))
    resEz = np.zeros((Np, K))
    resPEz = np.zeros((Np, K))

    # Cálculo do passo de tempo estável
    # JacobiGQ retorna os nós de quadratura r e os pesos w (não usados aqui)
    rLGL, _ = JacobiGQ(0, 0, N)
    # Distância mínima entre nós de Gauss-Lobatto
    rmin = abs(rLGL[1] - rLGL[0])
    # Ajusta dt de acordo com o fator de escala e rmin
    dtscale = dtscale3D()
    dt = 0.5 * np.min(dtscale) * rmin * 2.0 / 3.0

    #Para plotagem dos snapshots da solução
    t_snap = np.linspace(0, FinalTime, 6)
    isnap = 0
    #Inicializando as variáveis que vão salvar os dados para plotagem
    Ey_snaps = []
    x_snaps = []
    y_snaps = []
    t_record = []
    
    #Para monitorar e plotar o campo
    maskMonitor = np.where(
                np.isclose(x,0.2) &
                np.isclose(y,0.0) &
                np.isclose(z,0.2))
    nx, ny, nz = 21, 21, 21
    xi = np.linspace(-a, a, nx)
    yi = np.linspace(-b, b, ny)
    zi = np.linspace(-c, c, nz)
    Xi, Yi, Zi = np.meshgrid(xi, yi, zi)
    xMon = 0.2
    yMon = 0.2
    zMon = 0.3
    dist = np.sqrt((x.flatten('F')- xMon)**2 + (y.flatten('F')- yMon)**2 + (z.flatten('F')- zMon)**2)    
    i_closest = np.argmin(dist)
    Eplot = []
    EplotEx = []
    Hplot = []
    HplotEx = []
    time_vec = []



    # Laço principal de iteração de tempo
    while time < FinalTime:
        # Ajusta dt se estivermos próximos do tempo final
        if time + dt > FinalTime:
            dt = FinalTime - time

        # 5 estágios do Runge-Kutta de 4ª ordem
        for intrk in range(5):
            timelocal = time + rk4c[intrk] * dt

            # Calcula o lado direito das equações de Maxwell (modo TM)
            rhsHx, rhsHy, rhsHz, rhsEx, rhsEy, rhsEz, rhsPHx, rhsPHy, rhsPHz, rhsPEx, rhsPEy, rhsPEz = MaxwellRHS3D(Hx, Hy, Hz, Ex, Ey, Ez, PHx, PHy, PHz, PEx, PEy, PEz, timelocal)

            # Atualiza resíduos para cada componente
            resHx = rk4a[intrk] * resHx + dt * rhsHx
            resHy = rk4a[intrk] * resHy + dt * rhsHy
            resHz = rk4a[intrk] * resHz + dt * rhsHz
            resEx = rk4a[intrk] * resEx + dt * rhsEx
            resEy = rk4a[intrk] * resEy + dt * rhsEy
            resEz = rk4a[intrk] * resEz + dt * rhsEz 

            #Variáveis auxiliatórias da PML
            resPHx = rk4a[intrk] * resPHx + dt * rhsPHx
            resPHy = rk4a[intrk] * resPHy + dt * rhsPHy
            resPHz = rk4a[intrk] * resPHz + dt * rhsPHz
            resPEx = rk4a[intrk] * resPEx + dt * rhsPEx
            resPEy = rk4a[intrk] * resPEy + dt * rhsPEy
            resPEz = rk4a[intrk] * resPEz + dt * rhsPEz

            # Atualiza campos usando coeficientes de Runge-Kutta
            Hx = Hx + rk4b[intrk] * resHx
            Hy = Hy + rk4b[intrk] * resHy
            Hz = Hz + rk4b[intrk] * resHz
            Ex = Ex + rk4b[intrk] * resEx
            Ey = Ey + rk4b[intrk] * resEy
            Ez = Ez + rk4b[intrk] * resEz
            #Variáveis auxiliatórias ADE da PML
            PHx = PHx + rk4b[intrk] * resPHx
            PHy = PHy + rk4b[intrk] * resPHy
            PHz = PHz + rk4b[intrk] * resPHz
            PEx = PEx + rk4b[intrk] * resPEx
            PEy = PEy + rk4b[intrk] * resPEy
            PEz = PEz + rk4b[intrk] * resPEz

        #Salvando dados para plotar o campo Elétrico e o campo Magnético
        Eplot.append(Ey.flatten('F')[i_closest])
        Hplot.append(Hz.flatten('F')[i_closest])
        Eplotex = np.cos(omega * time - k * x.flatten('F')[i_closest])
        EplotEx.append(Eplotex)
        HplotEx.append(Eplotex) #O campo magnético tem a mesma expressão que o elétrico
        time_vec.append(time)
        

        # Incrementa o tempo
        time += dt

        #Para plotagem no plano z = 0
        if isnap < len(t_snap) and time >= t_snap[isnap]:
            #FAZER UMA GRID PARA INTERPOLAR COM X, Y E Z ANTES DE SELECIONAR PRA DEPOIS PLOTAR
            #Máscara que seleciona o plano desejado
            Ei = griddata(
                (x.flatten('F'), y.flatten('F'), z.flatten('F')),
                Ey.flatten('F'),
                (Xi, Yi, Zi),
                method='linear'
            )
            mask = np.abs(Zi - 0) < 1e-1 #Tolerância, ajustável
            
            #Extrai coordenadas e campo
            x_snaps.append(Xi[mask])
            y_snaps.append(Yi[mask])
            Ey_snaps.append(Ei[mask])
            t_record.append(time)
            isnap += 1
            #plot_solution(Xi[mask], Yi[mask], Ei[mask], f'Campo em z=0 em t= {time}')

    
    # --- pós-processamento: interpolação e plotagem ---
    #Calculando campos exatos
    #EplotEx = np.cos(omega * timelocal - k * x.flatten('F')[i_closest])
    Eyinc = np.cos(omega * timelocal - k * x.flatten('F'))
    Hzinc = np.cos(omega * timelocal - k * x.flatten('F'))

    # cria uma figura com 2 linhas x 1 coluna de subplots
    fig, axs = plt.subplots(2, 1, figsize=(8, 6), sharex=True)

    # primeiro subplot: E
    axs[0].plot(time_vec, Eplot,    label='Solução numérica')
    axs[0].plot(time_vec, EplotEx,  label='Solução analítica')
    axs[0].set_title(f'Campo elétrico em ({xMon}, {yMon}, {zMon})')
    axs[0].set_ylabel('$E_y$')
    axs[0].legend()
    axs[0].grid(True)

    # segundo subplot: H
    axs[1].plot(time_vec, Hplot,    label='Solução numérica')
    axs[1].plot(time_vec, HplotEx,  label='Solução analítica')
    axs[0].set_title(f'Campo magnético em ({xMon}, {yMon}, {zMon})')
    axs[1].set_xlabel('Tempo')
    axs[1].set_ylabel('$H_z$')
    axs[1].legend()
    axs[1].grid(True)

    # ajusta layout para não sobrepor títulos/eixos
    fig.tight_layout()

    # Cria figura 2x3
    vmin = Ey_snaps[0].min()
    vmax = Ey_snaps[0].max()
    fig, axes = plt.subplots(2, 3, figsize=(12, 8))
    for idx, ax in enumerate(axes.flat):
        # plot
        plot_solution_ax(x_snaps[idx], y_snaps[idx], Ey_snaps[idx], f't={t_record[idx]:.2f}', ax, vmin, vmax)
        #ax.set_aspect('equal')
        #ax.set_title(f"Ez @ z=0, t={t_record[idx]:.3f}")
        #fig.colorbar(cf, ax=ax, label=r'$E_z$')
    fig.tight_layout()
    return Hx, Hy, Hz, Ex, Ey, Ez, time

def plot_solution(x, y, u, label):

        x_flat = x.flatten('F')
        y_flat = y.flatten('F')
        u_flat = u.flatten('F')

        num_x, num_y = 200, 200
        xi = np.linspace(-a, a, num_x)
        yi = np.linspace(-b, b, num_y)
        Xi, Yi = np.meshgrid(xi, yi)

        # interpola Ez_flat nos pontos (x_flat, y_flat)
        Zi = griddata(
            points=(x_flat, y_flat),
            values=u_flat,
            xi=(Xi, Yi),
            method='cubic'
        )

        plt.figure(figsize=(6,5))
        cf = plt.contourf(Xi, Yi, Zi, levels=50)
        plt.colorbar(cf, label=r'$E_z$')
        plt.xlabel('x')
        plt.ylabel('y')
        plt.title(label)
        plt.axis('equal')

def plot_solution_ax(x, y, u, title, ax, vmin, vmax):
        x_flat = x.flatten('F')
        y_flat = y.flatten('F')
        u_flat = u.flatten('F')

        #grade regular para interpolação
        num_x, num_y = 200, 200
        xi = np.linspace(-a, a, num_x)
        yi = np.linspace(-b, b, num_y)
        Xi, Yi = np.meshgrid(xi, yi)

        #Interpolando
        Zi = griddata(
            points=(x_flat, y_flat),
            values=u_flat,
            xi = (Xi, Yi),
            method='cubic'
        )

        #faz o contourf no axes fornecido
        cf = ax.contourf(Xi, Yi, Zi, levels=50, vmin = vmin, vmax = vmax)
        ax.set_title(title)
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        #ax.set_xlim(-(a - Lpml - Ls), (a - Lpml - Ls))
        #ax.set_ylim(-(b - Lpml - Ls), (b - Lpml - Ls))
        ax.set_aspect('equal')
        #Adiciona colorbal individual
        plt.colorbar(cf, ax=ax, label=r'$E_z$')

        # Cria um retângulo para PML 
        rect = patches.Rectangle(
        (-a + Lpml, -b + Lpml),  # canto inferior esquerdo
        2 * (a - Lpml),              # largura
        2 * (b - Lpml),              # altura
        linewidth=1.5,
        edgecolor='black',
        facecolor='none',
        label = 'Domínio físico'
        )

        # Cria um retângulo para TFSF
        rect2 = patches.Rectangle(
        (-a + Lpml + Ls, -b + Lpml + Ls),  # canto inferior esquerdo
        2 * (a - Lpml - Ls),              # largura
        2 * (b - Lpml - Ls),              # altura
        linewidth=1.5,
        edgecolor='red',
        facecolor='none',
        label = 'Domínio TFSF'
        )
        ax.add_patch(rect2)
        ax.add_patch(rect)

def mesh_cubeK252(x0):
    h = 4
    L = 1

    gmsh.model.add("cubeK252_upml")
    factory = gmsh.model.occ

    # --- Criação da Geometria ---
    total_field_domain = factory.addBox(-x0, -x0, -x0, 2*x0, 2*x0, 2*x0)
    scattered_field_boundary = factory.addBox(-x0 - L, -x0 - L, -x0 - L, 2*(x0 + L), 2*(x0 + L), 2*(x0 + L))
    pml_boundary = factory.addBox(-x0 - 2*L, -x0 - 2*L, -x0 - 2*L, 2*(x0 + 2*L), 2*(x0 + 2*L), 2*(x0 + 2*L))

    pml_map, _ = factory.cut([(3, pml_boundary)], [(3, scattered_field_boundary)], removeTool=False)
    scattered_field_map, _ = factory.cut([(3, scattered_field_boundary)], [(3, total_field_domain)], removeTool=False)

    factory.synchronize()

    # --- Definição dos Grupos Físicos ---
    # Extrai a tag numérica de cada volume criado.
    # O resultado das operações de corte pode ter mais de um volume, mas aqui sabemos que é apenas um.
    pml_tag = pml_map[0][1]
    scattered_field_tag = scattered_field_map[0][1]
    # A tag do domínio total-field é a original, que geralmente é 1
    total_field_tag = total_field_domain

    # Atribui cada volume a um grupo físico com um nome e uma tag numérica única.
    # A dimensão é 3 para volumes.
    gmsh.model.addPhysicalGroup(3, [total_field_tag], tag=301, name="Total Field")
    gmsh.model.addPhysicalGroup(3, [scattered_field_tag], tag=302, name="Scattered Field")
    gmsh.model.addPhysicalGroup(3, [pml_tag], tag=303, name="PML")

    # --- Malha ---
    gmsh.option.setNumber("Mesh.MeshSizeMax", h)
    gmsh.option.setNumber("Mesh.MeshSizeMin", h)
    gmsh.model.mesh.generate(3)

def calculate_pml():
    #Matriz com as condutâncias do domínio
    sigma_x = np.zeros_like(x)
    sigma_y = np.zeros_like(y)
    sigma_z = np.zeros_like(z)
    #Em x
    maskXP = np.where(x > xpml)
    maskXN = np.where(x < -xpml) 
    sigma_x[maskXP] = sigma0x * ((x[maskXP]-xpml)**p)
    sigma_x[maskXN] = sigma0x * ((x[maskXN]+xpml)**p)
    #Em y
    maskYP = np.where(y > ypml)
    maskYN = np.where(y < -ypml)
    sigma_y[maskYP] = sigma0y * ((y[maskYP]-ypml)**p)
    sigma_y[maskYN] = sigma0y * ((y[maskYN]+ypml)**p)
    #Em z
    maskZP = np.where(z > zpml)
    maskZN = np.where(z < -zpml)
    sigma_z[maskZP] = sigma0z * ((z[maskZP]-zpml)**p)
    sigma_z[maskZN] = sigma0z * ((z[maskZN]+zpml)**p)

    print(f"Condutâncias PML:\nσx = {sigma_x}.")
    return sigma_x, sigma_y, sigma_z

def find_maskSource():
    # 1) decompõe vmapM em (linha, coluna)
    rows = vmapM % Np
    cols = vmapM // Np
    colsP =vmapP // Np

    # 2) pega coordenadas
    x_vals = x.flatten('F')[vmapM]
    y_vals = y.flatten('F')[vmapM]
    z_vals = z.flatten('F')[vmapM]

    # 3) máscara conjunta
    maskX = (
        (np.isclose(x_vals, xS) | np.isclose(x_vals, -xS))
        & (np.abs(y_vals) <= yS)
        & (np.abs(z_vals) <= zS)
    )

    # 4) contagem por coluna
    counts_per_colX = np.bincount(cols[maskX], minlength=K)

    # 5) colunas “válidas”
    valid_colsX = np.where(counts_per_colX > Nfp)[0]

    print("Colunas que têm ao menos", Nfp, 
        "pontos em x=xS com |y|≤yS e |z|≤zS:", valid_colsX)

    # 6) média de x em cada coluna (sobre todos os pontos do elemento)
    x_col_mean = x.mean(axis=0)   # vetor de tamanho K

    # 7) filtra colunas válidas de acordo com a média
    sct_cols_elemsX  = valid_colsX[np.abs(x_col_mean[valid_colsX]) > xS]
    total_cols_elemsX = valid_colsX[np.abs(x_col_mean[valid_colsX]) < xS]

    print("Colunas à esquerda (média x < -xS):",  sct_cols_elemsX)
    print("Colunas à direita (média x >  xS):",  total_cols_elemsX)

    #Para y
    # 3) máscara conjunta
    maskY = (
        (np.isclose(y_vals, yS) | np.isclose(y_vals, -yS))
        & (np.abs(x_vals) <= xS)
        & (np.abs(z_vals) <= zS)
    )

    # 4) contagem por coluna
    counts_per_colY = np.bincount(cols[maskY], minlength=K)

    # 5) colunas “válidas”
    valid_colsY = np.where(counts_per_colY > Nfp)[0]

    print("Colunas que têm ao menos", Nfp, 
        "pontos em y=yS com |x|≤xS e |z|≤zS:", valid_colsY)

    # 6) média de x em cada coluna (sobre todos os pontos do elemento)
    y_col_mean = y.mean(axis=0)   # vetor de tamanho K

    # 7) filtra colunas válidas de acordo com a média
    sct_cols_elemsY  = valid_colsY[np.abs(y_col_mean[valid_colsY]) > yS]
    total_cols_elemsY = valid_colsY[np.abs(y_col_mean[valid_colsY]) < yS]

    print("Colunas no lado espalhado (média |y| > yS):",  sct_cols_elemsY)
    print("Colunas no lado total (média |y| <  yS):",  total_cols_elemsY)

    #Para z
    # 3) máscara conjunta
    maskZ = (
        (np.isclose(z_vals, zS) | np.isclose(z_vals, -zS))
        & (np.abs(x_vals) <= xS)
        & (np.abs(y_vals) <= yS)
    )

    # 4) contagem por coluna
    counts_per_colZ = np.bincount(cols[maskZ], minlength=K)

    # 5) colunas “válidas”
    valid_colsZ = np.where(counts_per_colZ > Nfp)[0]

    print("Colunas que têm ao menos", Nfp, 
        "pontos em z=xS com |x|≤xS e |y|≤yS:", valid_colsZ)

    # 6) média de z em cada coluna (sobre todos os pontos do elemento)
    z_col_mean = z.mean(axis=0)   # vetor de tamanho K

    # 7) filtra colunas válidas de acordo com a média
    sct_cols_elemsZ  = valid_colsZ[np.abs(z_col_mean[valid_colsZ]) > zS]
    total_cols_elemsZ = valid_colsZ[np.abs(z_col_mean[valid_colsZ]) < zS]

    print("Colunas no lado espalhado (média |z| > zS):",  sct_cols_elemsZ)
    print("Colunas no lado total (média |z| <  zS):",  total_cols_elemsZ)


    maskScatterRegion = np.where((maskX & np.isin(cols, sct_cols_elemsX) & np.isin(colsP, total_cols_elemsX)
                                  | (maskY & np.isin(cols, sct_cols_elemsY) & np.isin(colsP, total_cols_elemsY)))
                                  | (maskZ & np.isin(cols, sct_cols_elemsZ) & np.isin(colsP, total_cols_elemsZ)))
    maskTotalregion = np.where((maskX & np.isin(cols, total_cols_elemsX) & np.isin(colsP, sct_cols_elemsX)
                                  | (maskY & np.isin(cols, total_cols_elemsY) & np.isin(colsP, sct_cols_elemsY)))
                                  | (maskZ & np.isin(cols, total_cols_elemsZ) & np.isin(colsP, sct_cols_elemsZ)))

    return maskScatterRegion, maskTotalregion

def Erro_l2(V, J, u, u_ex):

    """
    Computa a norma L2 do erro entre a solução numérica e exata
    em um esquema DG.

    Parâmetros
    ----------
    u : ndarray, forma (Np, K)
        Valores da solução numérica em cada nó/ponto de integração.
    V : ndarray, forma (Np, Np)
        Matriz de Vandermonde de referência.
    J : ndarray, forma (Np, K)
        Determinantes de Jacobiano (ou pesos de integração)
        para cada nó/ponto em cada elemento.

    Retorna
    -------
    errL2 : float
        Norma L2 global do erro.
    """
    # erro pontual em todos os elementos
    err = u_ex - u          # shape (Np, K)

    # matriz de massa de referência
    M = np.linalg.inv(V @ V.T)    # shape (Np, Np)

    # calcula quais elementos estão dentro do domínio físico
    # centro de massa de cada elemento
    x_mean = x.mean(axis=0)  # shape (K,)
    y_mean = y.mean(axis=0)
    z_mean = z.mean(axis=0)

    inside = (
        (np.abs(x_mean) < (a - Lpml - Ls)) &
        (np.abs(y_mean) < (b - Lpml - Ls)) &
        (np.abs(z_mean) < (c - Lpml - Ls))
    )                   # boolean mask shape (K,)

    # 4) para cada elemento dentro, monta o termo de erro local
    Np, K = err.shape
    err2_elem = np.zeros(K)
    for k in range(K):
        if not inside[k]:
            continue

        e  = err[:, k]       # vetor (Np,)
        Jk = J[:, k]         # vetor (Np,)
        # e^T * diag(Jk) * M * e = e • ( Jk * (M e) )
        err2_elem[k] = e @ (Jk * (M @ e))

    # 5) norma L2 global sobre os elementos filtrados
    errL2 = np.sqrt(np.sum(err2_elem))

    return errL2

N = 2 #Ordem polinomial
Np = (N+1) * (N+2) * (N+3) // 6 #Número de pontos
Nfp = (N+1) * (N+2) // 2
Nfaces = 4
NODETOL = 1e-10
alphaflux = 1

#Pontos no tetraedro equilátero
xeq, yeq, zeq = Nodes3D(N)

#Pontos no tetraedro de referência
r, s, t = xyztorst(xeq, yeq, zeq)

#Matrizes necessárias
V = Vandermonde3D(N, r, s, t)
Dr, Ds, Dt = DMatrices3D(N, r, s, t, V)
invV = np.linalg.inv(V)
MassMatrix = invV.T @ invV

#Criando a malha
#Parâmetros de dimensão: a, b, e c são os comprimentos em x, y e z a partir da origem

# gmsh.fltk.run()

lambda_ = 1 #Comprimento da onda
a = 3 #Tamanho de x a partir da origem (para + e -)
b = 3 #Tamanho de y a partir da origem (para + e -)
c = 3 #Tamanho de z a partir da origem (para + e -)
Lpml = lambda_ # Comprimento da PML
Ls = lambda_ #Comprimento da região TFSF
h = 16
refinement_level = np.sqrt(4 / np.sqrt(3)) / h
gmsh.initialize()
mesh_cubeK252(a - Lpml - Ls)
EToV = get_EToV(dim=3, index_based=0)
VX, VY, VZ = extract_VX_VY_VZ(get_nodes_data(dim=3))
# gmsh.fltk.run()
gmsh.finalize()
#vertices, EToV, edges, edge_map = generate_mesh(refinement_level, a, b, c, Lpml, Ls)
#export_mesh(vertices, EToV)
#VX = vertices[:, 0]
#VY = vertices[:, 1]
#VZ = vertices[:, 2]

np.save("C:\\git\\PyDG1D\\examplesData\\outputs\\cem_5\\fpb\\EToV.npy", EToV)
np.save("C:\\git\\PyDG1D\\examplesData\\outputs\\cem_5\\fpb\\VX.npy", VX)
np.save("C:\\git\\PyDG1D\\examplesData\\outputs\\cem_5\\fpb\\VY.npy", VY)
np.save("C:\\git\\PyDG1D\\examplesData\\outputs\\cem_5\\fpb\\VZ.npy", VZ)

Nv = len(VX)
K = len(EToV)
print(f'Número de tetraedros: {K}')
#plot_thetaedral_mesh(vertices, EToV, edges, edge_map)
#plot_tetrahedral_mesh2(VX, VY, VZ, EToV)

#Encontrando os pontos x, y e z na malha
va = EToV[:, 0]
vb = EToV[:, 1]
vc = EToV[:, 2]
vd = EToV[:, 3]
x = 0.5 * (-(1+r+s+t)[:, None] @ VX[None, va] + (1+r)[:, None] @ VX[None, vb] + (1+s)[:, None] @ VX[None, vc] + (1 + t)[:, None] @ VX[None, vd])
y = 0.5 * (-(1+r+s+t)[:, None] @ VY[None, va] + (1+r)[:, None] @ VY[None, vb] + (1+s)[:, None] @ VY[None, vc] + (1 + t)[:, None] @ VY[None, vd])
z = 0.5 * (-(1+r+s+t)[:, None] @ VZ[None, va] + (1+r)[:, None] @ VZ[None, vb] + (1+s)[:, None] @ VZ[None, vc] + (1 + t)[:, None] @ VZ[None, vd])

#Encontrando os pontos que estão nas interfaces (faces) dos tetraedros:
fmask1 = np.where(np.abs(t+1) < NODETOL)
fmask2 = np.where(np.abs(s+1) < NODETOL)
fmask3 = np.where(np.abs(r+s+t+1) < NODETOL)
fmask4 = np.where(np.abs(r+1) < NODETOL)
Fmask = np.concatenate([fmask1, fmask2, fmask3, fmask4]).T
Fx = x[Fmask.flatten('F')]
Fy = y[Fmask.flatten('F')]
Fz = z[Fmask.flatten('F')]

#Matriz de Lifting
LIFT = Lift3D(N, r, s, t)

#Fatores geométricos
rx, sx, tx, ry, sy, ty, rz, sz, tz, J = GeometricFactors3D(x, y, z, Dr, Ds, Dt)
nx, ny, nz, sJ = Normals3D()

#print(f'Jacobiano para verificação: {J}')
Fscale = sJ/J[Fmask.flatten('F')] #Por minha conta

#Matrizes e vetores de mapeamento
EToE, EToF = tiConnect3D(EToV)
mapM, mapP, vmapM, vmapP, vmapB, mapB = BuildMaps3D()

np.save("C:\\git\\PyDG1D\\examplesData\\outputs\\cem_5\\fpb\\mapM.npy", mapM)
np.save("C:\\git\\PyDG1D\\examplesData\\outputs\\cem_5\\fpb\\mapP.npy", mapP)
np.save("C:\\git\\PyDG1D\\examplesData\\outputs\\cem_5\\fpb\\vmapM.npy", vmapM)
np.save("C:\\git\\PyDG1D\\examplesData\\outputs\\cem_5\\fpb\\vmapP.npy", vmapP) 
np.save("C:\\git\\PyDG1D\\examplesData\\outputs\\cem_5\\fpb\\vmapB.npy", vmapB)
np.save("C:\\git\\PyDG1D\\examplesData\\outputs\\cem_5\\fpb\\mapB.npy", mapB)


print(f'vmapP: {vmapP}')
print(f'vmapM: {vmapM}')


#Computando operadores fracos
Vr, Vs, Vt = GradVandermonde3D(N, r, s, t)
Common_term = np.linalg.inv(V @ V.T)
Drw = np.linalg.solve((V @ V.T), (V @ Vr.T).T).T
Dsw = np.linalg.solve((V @ V.T), (V @ Vs.T).T).T
Dtw = np.linalg.solve((V @ V.T), (V @ Vt.T).T).T

#Coeficientes para Runge-Kutta
rk4a, rk4b, rk4c = rungekutta()

#Inicializando parâmetros da PML
xpml = a - Lpml
ypml = b - Lpml
zpml = c - Lpml
# sigma0x = 13.82
# sigma0y = 13.82
# sigma0z = 13.82
sigma0x = -np.log(1E-4)
sigma0y = -np.log(1E-4)
sigma0z = -np.log(1E-4)
p = 2
sigma_x, sigma_y, sigma_z = calculate_pml()

np.save("C:\\git\\PyDG1D\\examplesData\\outputs\\cem_5\\fpb\\sigma_x.npy", sigma_x)
np.save("C:\\git\\PyDG1D\\examplesData\\outputs\\cem_5\\fpb\\sigma_y.npy", sigma_y)
np.save("C:\\git\\PyDG1D\\examplesData\\outputs\\cem_5\\fpb\\sigma_z.npy", sigma_z)

#Encontrando as máscaras para fonte
threshold = N+2
xS = (a - Lpml - Ls)
yS = b - Lpml - Ls
zS = c - Lpml - Ls
maskScatterRegion, maskTotalRegion = find_maskSource()



#Resolvendo o problema
#Para inicializar o problema
#Parâmetros da fonte
k = 2 * np.pi / lambda_
f = 1 / lambda_
omega = 2 * np. pi * f
#Condição inicial
Hx = np.zeros((Np, K))
PHx = np.zeros((Np, K))
Hy = np.zeros((Np, K))
PHy = np.zeros((Np, K))
Hz = np.zeros((Np, K))
PHz = np.zeros((Np, K))
Ex = np.zeros((Np, K))
PEx = np.zeros((Np, K))
Ey = np.zeros((Np, K))
PEy = np.zeros((Np, K))
Ez = np.zeros((Np, K))
PEz = np.zeros((Np, K))

#Resolvendo o problema
FinalTime = 8

#Chama a rotina principal de integração temporal
Hx, Hy, Hz, Ex, Ey, Ez, time = Maxwell3D(Hx, Hy, Hz, Ex, Ey, Ez, PHx, PHy, PHz, PEx, PEy, PEz, FinalTime)

#Impressão de confirmação
print(f'Simulação concluída em t = {time:.3f}')

#Avaliando erro
#Para levantamento da norma L2 do erro
Sol_ex = np.cos(omega * time - k * x) #A solução exata é a mesma para o campo elétrico e magnético
#Definindo a máscara para o domínio físico
physical = (
        (np.abs(x) <= (a - Lpml - Ls)) &
        (np.abs(y) <= (b - Lpml - Ls)) &
        (np.abs(z) <= (c - Lpml - Ls))
    )                   # boolean mask shape (K,)
L2Ey = Erro_l2(V, J, Ey, Sol_ex)
print(f'O erro na norma L2 para o campo Ey foi de: {L2Ey}')
L2Hy = Erro_l2(V, J, Hz, Sol_ex)
print(f'O erro na norma L2 para o campo Hz foi de: {L2Hy}')

analit = np.cos(omega*time - k*x[physical])
msq = np.sqrt(1/(Np * K) * np.sum((Ey[physical] -analit)**2 ))
print(f'Erro quadratico médio: {msq}')
# plt.show()