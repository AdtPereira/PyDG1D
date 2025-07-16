import numpy as np
from scipy.linalg import eigh
from scipy.special import gamma
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import pandas as pd
from scipy.sparse import lil_matrix, eye
from math import ceil, pi
from mpl_toolkits.mplot3d import Axes3D
from meshpy.triangle import MeshInfo, build
from scipy.interpolate import griddata
from matplotlib.patches import Rectangle


import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as mtri

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

def JacobiGL(alpha, beta, Np):

    '''
    Função para computar os pontos e pesos de Legendgre-Gauss-Lobatto para no intervalo de [-1, 1]
    '''

    x = np.zeros((Np+1, 1))
    if Np==2:
        x[0] = -1
        x[1] = 1
        return x

    xint, w = JacobiGQ(alpha+1, beta+1, Np-2)
    x = np.concatenate(([-1], xint, [1]))
    return x

def Vandermonde1D(N, r):
    '''
    Função para calcular a matriz de Vandermonde:
    x - Pontos de entrada
    N - Grau do aproximação utilizado
    '''
    v1D = np.zeros((len(r), N + 1))
    for j in range(N + 1):
        v1D[:, j] = JacobiP(r[:], 0, 0, j)
    
    return v1D

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

def Warpfactor(N, rout):

    """
    Calcula a função warp de ordem N baseada nos nós de interpolação rout.

    Parameters
    ----------
    N : int
        Polynomial order.
    rout : array_like, shape (Nr,)
        Target points at which to evaluate the warp factor (in [-1, 1]).

    Returns
    -------
    warp : ndarray, shape (Nr,)
        The warp factor evaluated at each entry of rout.
    """
    # Compute LGL nodes (Jacobi–Gauss–Lobatto) and equidistant nodes
    LGLr = JacobiGL(0, 0, N + 1)                # shape: (N+1,)
    req = np.linspace(-1.0, 1.0, N + 1)     # shape: (N+1,)

    # Build the Vandermonde matrix at the equidistant nodes
    Veq = Vandermonde1D(N, req)             # shape: (N+1, N+1)

    # Number of target points
    Nr = len(rout)

    # Build Pmat by evaluating JacobiP at each rout for orders 0..N
    Pmat = np.zeros((N + 1, Nr))
    for i in range(N + 1):
        # In MATLAB: Pmat(i,:) = JacobiP(rout, 0, 0, i-1)'
        # In Python, i runs 0..N, so we pass k = i
        Pmat[i, :] = JacobiP(rout, 0, 0, i)

    # Solve for the Lagrange coefficient matrix Lmat:  Veq' \ Pmat
    # No MATLAB: Lmat = Veq' \ Pmat
    # Usando NumPy: solve (Veq.T) * X = Pmat  ⇒  X = inv(Veq.T) @ Pmat
    Lmat = np.linalg.solve(Veq.T, Pmat)     # shape: (N+1, Nr)

    # Compute the unscaled warp factor: warp = Lmat' * (LGLr - req)
    # (LGLr - req) has shape (N+1,), Lmat.T has shape (Nr, N+1)
    warp = Lmat.T.dot(LGLr - req)           # shape: (Nr,)

    # Aplicando o fator de escala. 
    # No MATLAB:
    #   zerof = (abs(rout) < 1.0 - 1.0e-10);
    #   sf    = 1.0 - (zerof.*rout).^2;
    #   warp  = warp./sf + warp.*(zerof-1);
    #
    # Fazendo em NumPy:
    zerof = np.abs(rout) < (1.0 - 1.0e-10)   # boolean mask, shape (Nr,)
    sf = 1.0 - (rout * zerof) ** 2           # shape: (Nr,)
    # Convert zerof to float for arithmetic: True→1.0, False→0.0
    zf_float = zerof.astype(float)
    warp = warp / sf + warp * (zf_float - 1.0)

    return warp

def Nodes2D(N):
    """
    Obtem nós (x, y) no triângulo equilátero para polinômios de ordem N.

    Parameters
    ----------
    N : int
        Polynomial order.
    
    Returns
    -------
    x : numpy.ndarray, shape (Np,)
        x-coordinates of the nodes.
    y : numpy.ndarray, shape (Np,)
        y-coordinates of the nodes.
    """
    # Parâmetros alpha otimizados para N<16
    alpopt = np.array([
        0.0000, 0.0000, 1.4152, 0.1001, 0.2751, 0.9800, 1.0999,
        1.2832, 1.3648, 1.4773, 1.4959, 1.5743, 1.5770, 1.6223, 1.6258
    ])
    
    # Configura os parâmetros alpha otimizados
    if N < 16:
        alpha = alpopt[N]
    else:
        alpha = 5.0 / 3.0
    
    # Número total de nós
    Np = int((N + 1) * (N + 2) / 2)

    # Aloca os arrays de coordenadas baricêntricas
    L1 = np.zeros(Np)
    L2 = np.zeros(Np)
    L3 = np.zeros(Np)
    
    # Preenche as coordenadas baricêntricas com a forma equidistante
    k = 0
    for n0 in range(N + 1):              
        for m0 in range(N + 1 - n0):     
            L1[k] = n0 / N
            L3[k] = m0 / N
            k += 1
    # Calcula L2 da relação L1 + L2 + L3 = 1
    L2 = 1.0 - L1 - L3
    
    # Converte as coordenadas baricêntricas para (x,y) no triângulo equilátero
    x = -L2 + L3
    y = (-L2 - L3 + 2.0 * L1) / np.sqrt(3.0)
    
    # Calcula as funções de blending para cada aresta
    blend1 = 4.0 * L2 * L3
    blend2 = 4.0 * L1 * L3
    blend3 = 4.0 * L1 * L2
    
    # Calcula os fatores de warp ao longo de cada aresta
    # warpf1 corresponde à aresta entre L3-L2, etc.
    warpf1 = Warpfactor(N, L3 - L2)
    warpf2 = Warpfactor(N, L1 - L3)
    warpf3 = Warpfactor(N, L2 - L1)
    
    # Combina blending & warp com o termo alpha
    warp1 = blend1 * warpf1 * (1.0 + (alpha * L1) ** 2)
    warp2 = blend2 * warpf2 * (1.0 + (alpha * L2) ** 2)
    warp3 = blend3 * warpf3 * (1.0 + (alpha * L3) ** 2)
    
    # Acumula as deformações associadas a cada aresta
    x += warp1 + np.cos(2.0 * np.pi / 3.0) * warp2 + np.cos(4.0 * np.pi / 3.0) * warp3
    y += 0.0 * warp1 + np.sin(2.0 * np.pi / 3.0) * warp2 + np.sin(4.0 * np.pi / 3.0) * warp3
    
    return x, y

def xytors(x, y):
    """
    Converte coordenadas (x, y) do triângulo equilátero para coordenadas (r, s) no
    triângulo de referência.

    Parameters
    ----------
    x : float or numpy.ndarray
        x-coordinate(s) in the equilateral triangle.
    y : float or numpy.ndarray
        y-coordinate(s) in the equilateral triangle.

    Returns
    -------
    r : float or numpy.ndarray
        r-coordinate(s) in the standard triangle.
    s : float or numpy.ndarray
        s-coordinate(s) in the standard triangle.
    """
    # Calcula as coordenadas baricêntricas L1, L2, L3 de (x, y)
    L1 = (np.sqrt(3.0) * y + 1.0) / 3.0
    L2 = (-3.0 * x - np.sqrt(3.0) * y + 2.0) / 6.0
    L3 = ( 3.0 * x - np.sqrt(3.0) * y + 2.0) / 6.0

    # Mapeia as coordenadas baricêntricas (L1, L2, L3) em (r, s) no triângulo de referência
    r = -L2 + L3 - L1
    s = -L2 - L3 + L1

    return r, s

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

def GradSimplex2DP(a, b, id, jd):
    """
    grad_simplex_2d_p: Retorna as derivadas da base modal (id, jd)
    no simplexo 2D nos pontos (a, b).

    Parâmetros:
    -----------
    a : float ou ndarray
        Coordenada a no domínio de referência.
    b : float ou ndarray
        Coordenada b no domínio de referência.
    id : int
        Grau do polinômio de Jacobi em a (índice modal).
    jd : int
        Grau do polinômio de Jacobi em b (índice modal).

    Retorna:
    --------
    dmodedr : float ou ndarray
        Derivada da base modal em relação a r (coordenada r).
    dmodeds : float ou ndarray
        Derivada da base modal em relação a s (coordenada s).
    """

    # Avalia o polinômio de Jacobi P_id^(0,0)(a) e sua derivada dP_id^(0,0)/da em a
    fa = JacobiP(a, 0, 0, id)
    dfa = GradJacobiP(a, 0, 0, id)

    # Avalia o polinômio de Jacobi P_jd^(2*id+1,0)(b) e sua derivada dP_jd^(2*id+1,0)/db em b
    gb = JacobiP(b, 2*id + 1, 0, jd)
    dgb = GradJacobiP(b, 2*id + 1, 0, jd)

    # ------------------------------------------------------
    # Cálculo da derivada em relação a r (d/dr)
    # Observação: r-derivative = (da/dr) * (d/da) + (db/dr) * (d/db)
    # Para transformações do simplexo, d/dr = (2/(1 - s)) * d/da = (2/(1 - b)) * d/da
    # ------------------------------------------------------
    # Primeiro termo: dfa * gb
    dmodedr = dfa * gb

    # Se id > 0, multiplica por ((0.5*(1 - b))^(id - 1))
    if id > 0:
        dmodedr = dmodedr * ((0.5 * (1 - b)) ** (id - 1))

    # Normalização final do termo r
    # Fator global: 2^(id + 0.5)
    dmodedr = (2 ** (id + 0.5)) * dmodedr

    # ------------------------------------------------------
    # Cálculo da derivada em relação a s (d/ds)
    # Observação: s-derivative = ((1 + a)/2)/((1 - b)/2) * d/da + d/db
    #             = ((1 + a)/(1 - b)) * d/da + d/db
    # ------------------------------------------------------
    # Primeiro termo: dfa * (gb * (0.5 * (1 + a)))
    dmodeds = dfa * (gb * (0.5 * (1 + a)))

    # Se id > 0, multiplica este termo parcial por ((0.5*(1 - b))^(id - 1))
    if id > 0:
        dmodeds = dmodeds * ((0.5 * (1 - b)) ** (id - 1))

    # Segundo termo: precisamos computar dgb * ((0.5*(1 - b))^id)
    tmp = dgb * ((0.5 * (1 - b)) ** id)

    # Se id > 0, subtrair o termo 0.5 * id * gb * ((0.5*(1 - b))^(id - 1))
    if id > 0:
        tmp = tmp - 0.5 * id * gb * ((0.5 * (1 - b)) ** (id - 1))

    # Soma o segundo termo multiplicado por fa
    dmodeds = dmodeds + fa * tmp

    # Normalização final do termo s
    dmodeds = (2 ** (id + 0.5)) * dmodeds

    return dmodedr, dmodeds

def GradVandermonde2D(N, r, s):
    """
    GradVandermonde2D: Inicializa o gradiente da base modal (i,j) em 2D nos pontos (r,s) para ordem N.
    
    Parâmetros de entrada:
    - N : inteiro
        A ordem máxima do polinômio modal.
    - r : array-like de formato (M,)
        Coordenadas r dos pontos de avaliação (M pontos no total).
    - s : array-like de formato (M,)
        Coordenadas s dos pontos de avaliação (M pontos no total).
    
    Retorna:
    - V2Dr : ndarray de forma (M, K)
        Matriz contendo as derivadas parciais em relação a r de cada função modal,
        onde K = (N+1)*(N+2)/2 é o número total de modos (i,j) com i+j ≤ N.
    - V2Ds : ndarray de forma (M, K)
        Matriz contendo as derivadas parciais em relação a s de cada função modal,
        mesmo padrão de índices que V2Dr.    
    """
    
    # Inicializa matrizes de derivadas com zeros
    V2Dr = np.zeros((len(r), int((N+1)*(N+2)/2)))
    V2Ds = np.zeros((len(r), int((N+1)*(N+2)/2)))
    
    # Converte (r,s) em coordenadas tensor-product (a,b) no simplex de referência
    # (supondo que rstoab esteja definido em outro módulo)
    a, b = rstoab(r, s)
    
    # Índice do modo corrente (coluna das matrizes V2Dr, V2Ds)
    sk = 0
    
    # Percorre todos os pares (i,j) tais que i+j ≤ N
    for i in range(0, N + 1):
        for j in range(0, N + 1 - i):
            # Chama a função que retorna as derivadas de uma função modal (i,j) nos pontos (a,b)
            # GradSimplex2DP deve devolver dois vetores de comprimento M:
            #   * dmod/dr em relação às coordenadas (r,s), mas internamente usa a relação (a,b)
            #   * dmod/ds em relação às coordenadas (r,s)
            dmod_r, dmod_s = GradSimplex2DP(a, b, i, j)
            
            # Armazena as derivadas nas colunas correspondentes
            V2Dr[:, sk] = dmod_r
            V2Ds[:, sk] = dmod_s
            
            # Avança para o próximo modo
            sk += 1
    
    return V2Dr, V2Ds

def Dmatrices2D(N, r, s, V):

    """
    Inicializa as matrizes de diferenciação (Dr, Ds) em coordenadas (r, s)
    no simplex, avaliadas nos pontos (r, s) para ordem N.

    Parâmetros:
    -----------
    N : int
        Ordem do polinômio (grau máximo).
    r : array_like, shape (Np,)
        Coordenadas dos pontos de interpolação na direção r.
    s : array_like, shape (Np,)
        Coordenadas dos pontos de interpolação na direção s.
    V : ndarray, shape (Np, Np)
        Matriz de Vandermonde 2D avaliada nos pontos (r, s) até ordem N.

    Retorna:
    --------
    Dr : ndarray, shape (Np, Np)
        Matriz de diferenciação em relação a r.
    Ds : ndarray, shape (Np, Np)
        Matriz de diferenciação em relação a s.
    """

    # Calcula as matrizes de gradiente da Vandermonde 2D em r e s
    Vr, Vs = GradVandermonde2D(N, r, s)

    # Em MATLAB: Dr = Vr / V  equivale a Vr * inv(V)
    # Em Python/Numpy, fazemos Dr = Vr.dot(inv(V)) ou usando solve para maior estabilidade
    # O mesmo vale para Ds
    V_inv = np.linalg.inv(V)  # inversa da matriz Vandermonde
    Dr = Vr.dot(V_inv)
    Ds = Vs.dot(V_inv)

    return Dr, Ds

def Filter2D(Norder, Nc, sp):

    """
    Inicializa a matriz de filtro 2D de ordem Norder, com corte em Nc e
    expoente de filtragem sp.
    
    Parâmetros:
    -----------
    Norder : int
        Ordem do polinômio (grau máximo) no espaço bidimensional.
    Nc : int
        Índice de corte (soma dos índices i+j a partir do qual o filtro é aplicado).
    sp : float ou int
        Expoente do filtro (grau da suavização exponencial).
    
    Retorna:
    --------
    F : ndarray de forma (Np, Np)
        Matriz de filtro 2D, onde Np = número de modos = (Norder+1)*(Norder+2)/2.
    """

    # Inicializa o vetor diagonal do filtro com 1.0 para todos os modos
    filterdiag = np.ones((Norder + 1) * (Norder + 2) // 2)

    # Constante alpha = -log(eps), onde eps é a precisão de ponto flutuante dupla
    alpha = -np.log(eps)  #Não entendi

    # Construção do filtro exponencial:
    # percorre todos os pares (i, j) tais que i+j <= Norder
    sk = 0  # índice linear no vetor filterdiag (0-based em Python)
    for i in range(Norder + 1):
        for j in range(Norder + 1 - i):
            # Se a soma dos índices de modo i+j for >= Nc, aplica o filtro
            if (i + j) >= Nc:
                # Fórmula do filtro exponencial:
                # exp(-alpha * [ ((i+j)-Nc)/(Norder-Nc) ]^sp )
                frac = ( (i + j) - Nc ) / (Norder - Nc)
                filterdiag[sk] = np.exp(-alpha * (frac ** sp))
            # Incrementa o índice linear
            sk += 1

    # Monta a matriz de filtro: F = V * diag(filterdiag) * invV
    # V e invV devem vir de Globals2D (já importados no topo)
    F = V @ np.diag(filterdiag) @ invV

    return F

def Lift2D():

    """
    Calcula o termo lift para o lado direito
    """

    # Initialize Emat as zero matrix of size (Np) × (Nfaces * Nfp)
    Emat = np.zeros((Np, Nfaces * Nfp))

    # ————— Face 1 —————
    # Select the reference‐r coords of the face‐1 nodes:
    faceR = r[Fmask[:, 0]]                    # shape: (Nfp,)
    V1D   = Vandermonde1D(N, faceR)           # shape: (Nfp, Np1D)
    # Build the 1D “mass‐edge” matrix: inv(V1D * V1Dᵀ)  → (Nfp × Nfp)
    massEdge1 = np.linalg.inv(V1D.dot(V1D.T))  # shape: (Nfp, Nfp)
    # Insert into Emat: rows = face‐1 node indices, cols = [0 : Nfp)
    Emat[Fmask[:, 0], 0 : Nfp] = massEdge1

    # ————— Face 2 —————
    faceR = r[Fmask[:, 1]]                    # shape: (Nfp,)
    V1D   = Vandermonde1D(N, faceR)           # shape: (Nfp, Np1D)
    massEdge2 = np.linalg.inv(V1D.dot(V1D.T))  # shape: (Nfp, Nfp)
    # Place into Emat at columns [Nfp : 2*Nfp)
    Emat[Fmask[:, 1], Nfp : 2 * Nfp] = massEdge2

    # ————— Face 3 —————
    # For the third face, we use the reference‐s coordinate:
    faceS = s[Fmask[:, 2]]                    # shape: (Nfp,)
    V1D   = Vandermonde1D(N, faceS)           # shape: (Nfp, Np1D)
    massEdge3 = np.linalg.inv(V1D.dot(V1D.T))  # shape: (Nfp, Nfp)
    Emat[Fmask[:, 2], 2 * Nfp : 3 * Nfp] = massEdge3

    # Finally, compute the “lift” operator:
    #    LIFT = V * (Vᵀ * Emat)
    LIFT = V.dot(V.T.dot(Emat))

    return LIFT

def Grad2D(u):
    """
    Função: Grad2D
    Propósito: Computar o campo gradiente 2D do escalar u.
    
    Parâmetro de entrada:
    u -- vetor coluna contendo os valores do escalar em cada nó (tamanho Np x 1)
    
    Saídas:
    ux -- componente x do gradiente (tamanho Npx1)
    uy -- componente y do gradiente (tamanho Npx1)
    """

    # Derivadas de u em coordenadas de referência (r, s)
    # Dr e Ds são matrizes de derivação nas direções r e s, respectivamente
    ur = Dr @ u       # du/dr
    us = Ds @ u       # du/ds

    # Mapeamento das derivadas de referência para coordenadas físicas (x, y)
    # rx, sx, ry, sy são arrays (de tamanho Np×1) contendo componentes do jacobiano inverso
    ux = rx * ur + sx * us    # du/dx = rx·(du/dr) + sx·(du/ds)
    uy = ry * ur + sy * us    # du/dy = ry·(du/dr) + sy·(du/ds)

    return ux, uy
    
def Div2D(u, v):

    """
    Calcula a divergência 2D do campo vetorial (u, v).
    
    Parâmetros de entrada:
    u -- componente x do campo vetorial (vetor coluna de tamanho Npx1)
    v -- componente y do campo vetorial (vetor coluna de tamanho Npx1)
    
    Saída:
    divu -- valor da divergência em cada nó (vetor coluna de tamanho Npx1)
    """

    # Derivadas parciais de u em coordenadas de referência (r, s)
    # Dr e Ds são matrizes de derivação nas direções r e s, respectivamente
    ur = Dr @ u     # du/dr
    us = Ds @ u     # du/ds

    # Derivadas parciais de v em coordenadas de referência (r, s)
    vr = Dr @ v     # dv/dr
    vs = Ds @ v     # dv/ds

    # Mapeamento das derivadas de referência para coordenadas físicas (x, y)
    # rx, sx, ry, sy são arrays (Np×1) contendo componentes do jacobiano inverso
    # A divergência em coordenadas físicas: du/dx + dv/dy
    divu = rx * ur + sx * us + ry * vr + sy * vs

    return divu

def Curl2D(ux, uy, uz):
    """
    Calcula o operador de rotacional 2D no plano (x, y).
    
    Parâmetros de entrada:
    ux -- componente x do campo vetorial (vetor coluna de tamanho Npx1)
    uy -- componente y do campo vetorial (vetor coluna de tamanho Npx1)
    uz -- componente z do campo vetorial (vetor coluna de tamanho Npx1). 
          Se passar vetor vazio, só será calculada a componente z do rotacional.
    
    Saídas:
    vx -- componente x do rotacional (vetor coluna de tamanho Npx1). Vazio se uz for vazio.
    vy -- componente y do rotacional (vetor coluna de tamanho Npx1). Vazio se uz for vazio.
    vz -- componente z do rotacional (vetor coluna de tamanho Npx1)
    """

    # Derivadas de ux em coordenadas de referência (r, s)
    uxr = Dr @ ux    # dux/dr
    uxs = Ds @ ux    # dux/ds

    # Derivadas de uy em coordenadas de referência (r, s)
    uyr = Dr @ uy    # duy/dr
    uys = Ds @ uy    # duy/ds

    # Cálculo da componente z do rotacional (∂uy/∂x - ∂ux/∂y)
    # Mapeamento de derivadas de referência para físicas:
    # vz = rx·(∂uy/∂r) + sx·(∂uy/∂s) - [ry·(∂ux/∂r) + sy·(∂ux/∂s)]
    vz = rx * uyr + sx * uys - ry * uxr - sy * uxs

    # Inicializa vx e vy como vazios
    vx = np.array([])  
    vy = np.array([])

    # Se uz não for vetor vazio, calcula componentes x e y do rotacional
    if uz.size != 0:
        # Derivadas de uz em coordenadas de referência (r, s)
        uzr = Dr @ uz   # ∂uz/∂r
        uzs = Ds @ uz   # ∂uz/∂s

        # Cálculo das componentes x e y do rotacional:
        # vx =  ry·(∂uz/∂r) + sy·(∂uz/∂s)   (∂uz/∂y)
        # vy = -rx·(∂uz/∂r) - sx·(∂uz/∂s)   (−∂uz/∂x)
        vx = ry * uzr + sy * uzs
        vy = -rx * uzr - sx * uzs

    return vx, vy, vz

def generate_mesh(target_elements, a, b, x0, y0, xpml, ypml, xMeta):

    mesh_info = MeshInfo()  # Instância de MeshInfo para criação da malha

    # 1) Pontos externos
    outer = [(-a, -b),  # 0
             ( a, -b),  # 1
             ( a,  b),  # 2
             (-a,  b)]  # 3
    # 2) Pontos médios da interface entre a PML e a TFSF region
    '''medium = [(-xpml, -ypml), #4
              (xpml, -ypml),  #5
              (xpml, ypml),   #6
              (-xpml, ypml)]  #7'''
    # 2) Pontos internos, centralizados em (0,0)
    ix, iy = x0 / 2.0, y0 / 2.0
    inner = [( ix, -iy),  # 4
             ( ix,  iy),  # 5
             (-ix,  iy),  # 6
             (-ix, -iy)]  # 7
    # 3) Pontos da metasuperfície
    meta = [(xMeta, -b), #8
            (xMeta, b)]  #9
    

    # Concatena todos os pontos
    #points = outer + medium + inner + meta
    points = outer + inner + meta
    mesh_info.set_points(points)

    # Define as facets (arestas) das fronteiras externa e interna
    # Facets externas
    facets = [ [0, 1], [1, 2], [2, 3], [3, 0] ]
    # Facets médias
    facets += [ [4, 5], [5, 6], [6, 7], [7, 4] ]
    '''# Facets internas
    facets += [ [8, 9], [9, 10], [10, 11], [11, 8] ]'''
    #Facets da metasuperfície
    facets += [ [12, 13], [13, 12]]
    #Facets para inserir a fonte
    mesh_info.set_facets(facets)

    # Calcula área total do domínio externo (para estimar max_area)
    total_area = 4 * a * b
    max_vol = total_area / float(target_elements)

    # Constrói a malha, preservando as facets internas e externas
    mesh = build(mesh_info, max_volume=max_vol)

    # Extrai vértices e elementos
    vertices = np.array(mesh.points)
    elements = np.array(mesh.elements)

    return vertices, elements

def plot_mesh(vertices, elements, order):

    plt.figure(figsize=(8, 8))
    
    # Plotar os elementos
    for elem in elements:
        elem_vertices = vertices[elem[:3], :]  # Pegando os vértices principais do triângulo
        plt.fill(elem_vertices[:, 0], elem_vertices[:, 1], edgecolor='black', fill=False, linewidth=0.5)

        # Conexões entre os nós adicionais, se ordem > 1
        if order > 1:
            for i, vi in enumerate(elem):
                for vj in elem[i+1:]:
                    plt.plot(
                        [vertices[vi, 0], vertices[vj, 0]], 
                        [vertices[vi, 1], vertices[vj, 1]], 
                        color='gray', linestyle='dotted', linewidth=0.5
                    )

    # Plotar os nós
    plt.scatter(vertices[:, 0], vertices[:, 1], c='red', s=10, label='Nós')
    
    # Anotar os índices dos nós
    for i, (x, y) in enumerate(vertices):
        plt.text(x, y, f'{i}', fontsize=8, ha='right', va='top')

    plt.gca().set_aspect('equal', adjustable='box')
    plt.title(f'Malha de Elementos Finitos - Ordem {order}')
    plt.xlabel('r')
    plt.ylabel('z')
    plt.legend()
    plt.grid(True, linestyle='--', alpha=0.7)

def GeometricFactors2D(x, y, Dr, Ds):
    """
    Calcula os elementos métricos para o mapeamento local dos elementos.
    
    Parâmetros de entrada:
    x -- coordenadas x dos pontos da grid (Np x K)
    y -- coordenadas y dos pontos da grid (Np x K)
    Dr -- Matriz de derivadas em r
    Ds -- Matriz de derivadas em s
    
    Saídas:
    rx -- Fator de escala x em relação a r.
    sx -- Fator de escala x em relação a s.
    ry -- Fator de escala y em relação a r.
    sy -- Fator de escala y em relação a s.
    J  -- Jacobiano da transformação
    """

    # Calcula as derivadas de x em relação a r e s
    xr = Dr @ x
    xs = Ds @ x

    # Calcula as derivadas de y em relação a r e s
    yr = Dr @ y
    ys = Ds @ y

    # Determinante do Jacobiano: J = ∂(x,y)/∂(r,s) = xr*ys - xs*yr (com sinal trocado)
    J = -xs * yr + xr*ys

    # Fatores métricos: 
    rx = ys / J
    sx = -yr / J
    ry = -xs / J
    sy = xr/J

    return rx, sx, ry, sy, J

def Normals2D():

    """
    Calcula as normais apontando para fora nos faces dos elementos e os jacobianos de superfície.
    Retorna:
      nx (np.ndarray): componentes x das normais nos nós das faces (dimensão 3*Nfp x K)
      ny (np.ndarray): componentes y das normais nos nós das faces (dimensão 3*Nfp x K)
      sJ (np.ndarray): jacobiano de superfície (comprimento das normais não normalizadas)
    """

    #Calculando os vatores geométricos no interior de cada elemento
    xr = Dr @ x
    yr = Dr @ y
    xs = Ds @ x
    ys = Ds @ y
    Fmask2 = Fmask.flatten('F')

    #Interpolando os fatores geométricos para os nós das faces
    fxr = xr[Fmask2, :]
    fxs = xs[Fmask2, :]
    fyr = yr[Fmask2, :]
    fys = ys[Fmask2, :]

    #Inicializando as matrizes de normais não normalizadas para cada face
    nx = np.zeros((3 * Nfp, K))
    ny = np.zeros((3 * Nfp, K))

    #Definindo o índice para cada bloco de Nfp linhas:
    fid1 = np.linspace(0, Nfp-1, Nfp, dtype=int)
    fid2 = fid1 + Nfp
    fid3 = fid2 + Nfp

    #Calculando as componentes normais não normalizadas em cada face
    # Face 1
    nx[fid1, :] = fyr[fid1, :]      # componente x = dy/dr
    ny[fid1, :] = -fxr[fid1, :]     # componente y = -dx/dr

    # Face 2
    nx[fid2, :] = fys[fid2, :] - fyr[fid2, :]   # componente x = dy/ds - dy/dr
    ny[fid2, :] = -fxs[fid2, :] + fxr[fid2, :]   # componente y = -dx/ds + dx/dr

    # Face 3
    nx[fid3, :] = -fys[fid3, :]    # componente x = -dy/ds
    ny[fid3, :] = fxs[fid3, :]     # componente y = dx/ds

    # Normalizando as normais e calculando o jacobiano de superfície sJ, que é 
    # o comprimento da normal não-normalizada.
    sJ = np.sqrt(nx ** 2 + ny ** 2)   # comprimento da normal em cada ponto de face
    nx = nx / sJ                       # normalização: nx / |n|
    ny = ny / sJ                       # normalização: ny / |n|

    return nx, ny, sJ

def Connect2D(EToV):
    """
    Constroi a conectividade global para grid baseada em EToV padrão do gerador de grid.
    Entrada:
        EToV: array inteiro de forma (K, 3) com os índices dos nós de cada elemento triangular.
    Saída:
        EToE: array inteiro de forma (K, 3), onde EToE[k,f] é o elemento vizinho ao elemento k na face f.
        EToF: array inteiro de forma (K, 3), onde EToF[k,f] é o número da face de EToE[k,f] que faz fronteira com k.
    Tudo documentado em Português.
    """
    # Número de faces por elemento (triângulo)
    Nfaces = 3

    # Número de elementos e número máximo de vértices
    K = EToV.shape[0]
    Nv = int(EToV.max()) + 1

    # Total de faces na malha (cada elemento tem 3 faces)
    TotalFaces = Nfaces * K

    # Mapeamento local: para cada face local, quais são os (dois) vértices locais do elemento
    vn = np.array([[0, 1],
                   [1, 2],
                   [0, 2]], dtype=int)

    # Construir a matriz esparsa global "face -> vértice"
    # Linhas: cada face global (1..TotalFaces)
    # Colunas: cada vértice (1..Nv)
    # Em MATLAB, usou-se spalloc(TotalFaces, Nv, 2*TotalFaces) porque cada face se conecta a exatamente 2 vértices.
    SpFToV = lil_matrix((TotalFaces, Nv), dtype=int)

    # Percorre cada elemento e cada face local, preenchendo SpFToV
    sk = 0
    for k in range(K):
        for face in range(Nfaces):
            # Encontrar índices globais dos 2 nós que definem esta face local
            # EToV[k, vn[face,:]] no MATLAB; em Python, subtrai-se 1 se EToV for 1-base (supomos 1-base vinda do gerador)
            v1, v2 = EToV[k, vn[face, 0]], EToV[k, vn[face, 1]]
            SpFToV[sk, v1] = 1
            SpFToV[sk, v2] = 1
            sk += 1

    # Converter para formato CSR para multiplicação rápida
    SpFToV = SpFToV.tocsr()

    # Montar a matriz "face -> face" (esparsa)
    # SpFToF = SpFToV * SpFToV' - 2*I
    # Se duas faces compartilham exatamente os mesmos 2 nós, o produto terá valor 2 nessa posição.
    SpFToF = SpFToV.dot(SpFToV.transpose()) - 2 * eye(TotalFaces, dtype=int)

    # Encontrar pares de faces cujo valor em SpFToF é igual a 2 (significa que compartilham exatamente 2 nós)
    SpFToF_coo = SpFToF.tocoo()
    mask = (SpFToF_coo.data == 2)
    faces1 = SpFToF_coo.row[mask]  # índice da primeira face global
    faces2 = SpFToF_coo.col[mask]  # índice da segunda face global

    # Converter o índice global de face (1..TotalFaces) em par (elemento, face local)
    # Em MATLAB: element1 = floor((faces1-1)/Nfaces) + 1; face1 = mod((faces1-1),Nfaces) + 1;
    element1 = np.floor_divide(faces1, Nfaces)        # já zero-based
    face1    = np.mod(faces1, Nfaces)                # zero-based
    element2 = np.floor_divide(faces2, Nfaces)
    face2    = np.mod(faces2, Nfaces)

    # Inicializar EToE e EToF com valores padrão:
    # EToE[k,f] = k (vizinho consigo mesmo) e EToF[k,f] = f (mesma face), para depois sobrescrever
    EToE = np.tile(np.arange(K, dtype=int).reshape((K, 1)), (1, Nfaces))
    EToF = np.tile(np.arange(Nfaces, dtype=int).reshape((1, Nfaces)), (K, 1))

    # Índices de inserção: transforma (element, face) em índice único no array achatado
    # Em MATLAB usaria sub2ind([K, Nfaces], element1, face1), mas aqui fazemos manualmente
    ind = element1 * Nfaces + face1  # corresponde à posição linear na matriz K x Nfaces

    # Substituir nos locais apropriados: para cada par de faces encontradas
    EToE_flat = EToE.flatten()
    EToF_flat = EToF.flatten()
    EToE_flat[ind] = element2
    EToF_flat[ind] = face2

    # Voltar a redimensionar para (K, Nfaces)
    EToE = EToE_flat.reshape((K, Nfaces))
    EToF = EToF_flat.reshape((K, Nfaces))

    return EToE, EToF

def BuildMaps2D():

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

            # Cálculo da distância de referência entre vértices da face (usada na tolerância)
            v1 = EToV[k1, f1]
            v2 = EToV[k1, (f1 + 1) % Nfaces]
            refd = np.sqrt((VX[v1] - VX[v2])**2 + (VY[v1] - VY[v2])**2)

            # Índices globais dos nós da face atual (menor lado)
            vidM = vmapM[:, f1, k1]
            # Índices globais dos nós da face do vizinho (maior lado)
            vidP = vmapM[:, f2, k2]

            # Coordenadas dos nós das duas faces
            x1 = x.flatten('F')[vidM]
            y1 = y.flatten('F')[vidM]
            x2 = x.flatten('F')[vidP]
            y2 = y.flatten('F')[vidP]
            x1, y1 = x1[:, None] @ one, y1[:, None] @ one
            x2, y2 = x2[:, None] @ one, y2[:, None] @ one

            # Matriz de distâncias quadradas entre nós das duas faces
            D = (x1 - x2.T)**2 + (y1 - y2.T)**2

            # Identifica pares de nós que coincidem (dentro da tolerância)
            idM, idP = np.where(np.sqrt(np.abs(D)) < NODETOL * refd)

            # Atualiza os mapeamentos com os índices corretos dos vizinhos
            vmapP[idM, f1, k1] = vidP[idP]
            mapP[idM, f1, k1] = idP + f2 * Nfp + k2 * Nfaces * Nfp

    # Achata os arrays tridimensionais para vetores coluna
    vmapM = vmapM.flatten(order='F')
    vmapP = vmapP.flatten(order='F')
    mapM = mapM.flatten(order='F')
    mapP = mapP.flatten(order='F')

    # Identifica os nós de fronteira: onde o vizinho é o próprio nó (i.e., vmapP == vmapM)
    mapB = np.where(vmapP == vmapM)[0]
    vmapB = vmapM[mapB]

    return mapM, mapP, vmapM, vmapP, vmapB, mapB

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

def dtscale2D():
    """
    Faz o escalonamento temporal para fins de estabilidade
    """

    vmask1 = np.where(abs(s + r + 2) < NODETOL)[0]
    vmask2 = np.where(abs(r - 1) < NODETOL)[0]
    vmask3 = np.where(abs(s -1) < NODETOL)[0]
    vmask = np.concatenate([vmask1, vmask2, vmask3])
    vx, vy = x[vmask], y[vmask]

    #Calculando semi-perímetro e área
    len1 = np.sqrt((vx[0] - vx[1])**2 + (vy[0] - vy[1])**2)
    len2 = np.sqrt((vx[1] - vx[2])**2 + (vy[1] - vy[2])**2)
    len3 = np.sqrt((vx[2] - vx[0])**2 + (vy[2] - vy[0])**2)
    sper = (len1 + len2 + len3)/2
    Area = np.sqrt(sper * (sper-len1) * (sper-len2) * (sper-len3))
    #Calcula a escala usando o raio do círculo inscrito
    dtscale = Area/sper

    return dtscale

'''
FUNÇÃO PARA TRABALHAR AS CONDIÇÕES DE CONTORNO EM DIFERENTES NÓS
def BuildBCMaps2D(BCType):

    """
    Esta função constrói mapeamento nodal especializado para os diferentes tipos de condições
     de contorno, especificados pela variável BCType
     
    BCType deve ser uma matriz (Kx3) correspondente a uma entrada para cada face na grid. Aquelas que tiverem uma condição de contorno
    específica deverão ter valores diferentes de 0.

    Essa informação deve ser dada pelo gerador da grid

    No exemplo, considerou-se tipos diferentes de contorno com a designação de, por exemplo, "In" e "Out"
        In está associado ao valor 1
        Out está associado ao valor 2
    """

    bct = BCType.T
    bnodes = np.ones((Nfp, 1)) @ (bct.flatten('F')[None, :])
    bnodes = bnodes.flatten("F")
    maps, vmaps = {}, {}

    for i, BC in enumerate(BCs.keys()):
        maps[BC] = np.where(bnodes == BCs[BC])[0]
        vmaps[BC] = vmapM[maps[BC]]

'''

def MaxwellRHS2D(Hx, Hy, Ez, Px, Py, Pz, timelocal):
    """
    Avalia o lado direito da equação de Maxwell no caso TM.

    Parameters:
    -----------
    Hx, Hy, Ez : ndarray
        Field arrays of shape (Np, K)
    alpha : float, optional
        Upwind parameter (default is 1.0)

    Returns:
    --------
    rhsHx, rhsHy, rhsEz : ndarray
        Right-hand side arrays of shape (Np, K)
    """
    # Field differences at faces
    #dHx = np.zeros((Nfp * Nfaces, K))
    #dHy = np.zeros((Nfp * Nfaces, K))
    #dEz = np.zeros((Nfp * Nfaces, K))

    Hxflat = Hx.flatten('F')
    Hyflat = Hy.flatten('F')
    Ezflat = Ez.flatten('F')

    dHx = Hxflat[vmapM] - Hxflat[vmapP]
    dHy = Hyflat[vmapM] - Hyflat[vmapP]
    dEz = Ezflat[vmapM] - Ezflat[vmapP]

    if (contour == 1) :
        # Reflective boundary conditions (Ez+ = -Ez-) (Caso PEC)
        dHx[mapB] = 0.0
        dHy[mapB] = 0.0
        dEz[mapB] = 2.0 * Ezflat[vmapB]

    nxflat = nx.flatten('F')
    nyflat = ny.flatten('F')

    if metaOn == 1:
        #Para a metasuperfície, aplicando as condições de contorno generalizadas
        Ezi = Ezflat[vmapM[maskMetaleft]]
        Ezj = Ezflat[vmapM[maskMetaright]]
        Hyi = Hyflat[vmapM[maskMetaleft]]
        Hyj = Hyflat[vmapM[maskMetaright]] 
        #Coeficientes
        Xmm = - 2* (T - 1 - R) / (1 + T - R)
        Xee = - 2* (T - 1 - R) / (1 + T + R)
        #Calculando E0-, E0+, H0- e H0+
        E0neg = Ezj - Xmm * (Hyi + Hyj) / 2
        E0pos = Ezi + Xmm * (Hyi + Hyj) / 2
        H0neg = Hyj - Xee * (Ezi + Ezj) / 2
        H0pos = Hyi + Xee * (Ezi + Ezj) / 2
        #Uma outra opção seria usar vmapP para calcular os E0-, E0+, H0-, H0+ para combinar caso ache um número diferente de termos dos dois lados (maskMetaleft e maskMetaright)
        #Calculando Delta(Ez) e Delta(Hy) nas arestas virtuais
        dEz[maskMetaleft] = Ezflat[vmapM[maskMetaleft]] - E0neg
        dEz[maskMetaright] = Ezflat[vmapM[maskMetaright]] - E0pos
        dHy[maskMetaleft] = Hyflat[vmapM[maskMetaleft]] - H0neg
        dHy[maskMetaright] = Hyflat[vmapM[maskMetaright]] - H0pos

    #Agora para a fonte usando formulação TFSF
    #Calculando o campo incidente
    if (timelocal < FinalTime) & (SourceTFSF == 1):
        Ezinc = +E0 * np.sin(omega * timelocal - k * x.flatten('F'))
        Hyinc = -E0 * np.sin(omega * timelocal - k * x.flatten('F'))
        dEz[maskScatterRegion] = Ezflat[vmapM[maskScatterRegion]] - (Ezflat[vmapP[maskScatterRegion]] - Ezinc[vmapP[maskScatterRegion]])
        dHy[maskScatterRegion] = Hyflat[vmapM[maskScatterRegion]] - (Hyflat[vmapP[maskScatterRegion]] - Hyinc[vmapP[maskScatterRegion]])            
        
        dEz[maskTotalRegion] = Ezflat[vmapM[maskTotalRegion]] - (Ezflat[vmapP[maskTotalRegion]] + Ezinc[vmapP[maskTotalRegion]])
        dHy[maskTotalRegion] = Hyflat[vmapM[maskTotalRegion]] - (Hyflat[vmapP[maskTotalRegion]] + Hyinc[vmapP[maskTotalRegion]])
 
    # Upwind fluxes
    ndotdH = nxflat * dHx + nyflat * dHy
    fluxHx = nyflat * dEz + alpha * (ndotdH * nxflat - dHx)
    fluxHy = -nxflat * dEz + alpha * (ndotdH * nyflat - dHy)
    fluxEz = -nxflat * dHy + nyflat * dHx - alpha * dEz

    #Colocando no formato correto
    fluxHx = fluxHx.reshape((Nfp * Nfaces, K) , order='F')
    fluxHy = fluxHy.reshape((Nfp * Nfaces, K) , order='F')
    fluxEz = fluxEz.reshape((Nfp * Nfaces, K) , order='F')


    # Local derivatives of fields
    Ezx, Ezy = Grad2D(Ez)
    Hz = np.zeros_like(Hx)
    CuHx, CuHy, CuHz = Curl2D(Hx, Hy, Hz)

    # Compute right-hand sides
    rhsHx = -Ezy + LIFT @ (Fscale * fluxHx) / 2.0 + Px - (sigma_y - sigma_x) * Hx
    rhsHy = Ezx  + LIFT @ (Fscale * fluxHy) / 2.0 + Py - (sigma_x - sigma_y) * Hy
    rhsEz = CuHz + LIFT @ (Fscale * fluxEz) / 2.0 - Pz - (sigma_x + sigma_y) * Ez
    rhsPx = -sigma_x * Px + (sigma_x * sigma_y - sigma_x**2) * Hx
    rhsPy = -sigma_y * Py + (sigma_x * sigma_y - sigma_y**2) * Hy
    rhsPz = sigma_x * sigma_y * Ez

    return rhsHx, rhsHy, rhsEz, rhsPx, rhsPy, rhsPz

def Maxwell2D(Hx, Hy, Ez, Px, Py, Pz, FinalTime):
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
    resHy = np.zeros((Np, K))
    resEz = np.zeros((Np, K))
    
    #Resíduos de Px, Py, Qx, Qy
    resPx = np.zeros((Np, K))
    resPy = np.zeros((Np, K))
    resPz = np.zeros((Np, K))

    # Cálculo do passo de tempo estável
    # JacobiGQ retorna os nós de quadratura r e os pesos w (não usados aqui)
    rLGL, _ = JacobiGQ(0, 0, N)
    # Distância mínima entre nós de Gauss-Lobatto
    rmin = abs(rLGL[1] - rLGL[0])
    # Ajusta dt de acordo com o fator de escala e rmin
    dtscale = dtscale2D()
    dt = 0.05 * np.min(dtscale) * rmin * 2.0 / 3.0
    print(f'dt calculado: {dt}')

    #Para plotagem dos snapshots da solução
    t_snap = np.linspace(0, FinalTime, 6)
    isnap = 0
    #Inicializando as variáveis que vão salvar os dados para plotagem
    Ey_snaps = []
    x_snaps = []
    y_snaps = []
    t_record = []

    #Para monitorar e plotar o campo
    nx, ny = 21, 21
    xi = np.linspace(-a, a, nx)
    yi = np.linspace(-b, b, ny)
    Xi, Yi = np.meshgrid(xi, yi)
    xMon = 0.1
    yMon = 0
    xMon2 = -1
    dist = np.sqrt((x.flatten('F')- xMon)**2 + (y.flatten('F')- yMon)**2)   
    dist2 = np.sqrt((x.flatten('F')- xMon2)**2 + (y.flatten('F')- yMon)**2)    
    i_closest = np.argmin(dist)
    i_closest2 = np.argmin(dist2)
    Eplot = []
    Eplot2 = []
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
            rhsHx, rhsHy, rhsEz, rhsPx, rhsPy, rhsPz = MaxwellRHS2D(Hx, Hy, Ez, Px, Py, Pz, timelocal)

            # Atualiza resíduos para cada componente
            resHx = rk4a[intrk] * resHx + dt * rhsHx
            resHy = rk4a[intrk] * resHy + dt * rhsHy
            resEz = rk4a[intrk] * resEz + dt * rhsEz
            
            #Variáveis auxiliatórias da PML
            resPx = rk4a[intrk] * resPx + dt * rhsPx
            resPy = rk4a[intrk] * resPy + dt * rhsPy
            resPz = rk4a[intrk] * resPz + dt * rhsPz

            # Atualiza campos usando coeficientes de Runge-Kutta
            Hx = Hx + rk4b[intrk] * resHx
            Hy = Hy + rk4b[intrk] * resHy
            Ez = Ez + rk4b[intrk] * resEz
            
            #Variáveis auxiliatórias da PML
            Px = Px + rk4b[intrk] * resPx
            Py = Py + rk4b[intrk] * resPy
            Pz = Pz + rk4b[intrk] * resPz

        #Salvando dados para plotar o campo Elétrico e o campo Magnético
        Eplot.append(Ez.flatten('F')[i_closest])
        Hplot.append(Hy.flatten('F')[i_closest])
        Eplot2.append(Ez.flatten('F')[i_closest2])
        #Eplotex = np.sin(omega * timelocal - k * x.flatten('F')[i_closest])
        #EplotEx.append(Eplotex)
        #HplotEx.append(-Eplotex) #O campo magnético tem a mesma expressão que o elétrico
        time_vec.append(time)

        # Incrementa o tempo
        time += dt

        #Para plotagem dos snaps
        if isnap < len(t_snap) and time >= t_snap[isnap]:
            #FAZER UMA GRID PARA INTERPOLAR COM X, Y E Z ANTES DE SELECIONAR PRA DEPOIS PLOTAR
            #Máscara que seleciona o plano desejado
            Ei = griddata(
                (x.flatten('F'), y.flatten('F')),
                Ez.flatten('F'),
                (Xi, Yi),
                method='linear'
            )
            
            #Extrai coordenadas e campo
            x_snaps.append(Xi)
            y_snaps.append(Yi)
            Ey_snaps.append(Ei)
            t_record.append(time)
            isnap += 1

    # cria uma figura com 2 linhas x 1 coluna de subplots
    fig, axs = plt.subplots(2, 1, figsize=(8, 6), sharex=True)

    # primeiro subplot: E
    axs[0].plot(time_vec, Eplot2,    label='Solução numérica')
    #axs[0].plot(time_vec, EplotEx,  label='Solução analítica')
    axs[0].set_title(f'Electrical field in ({xMon2}, {yMon})')
    axs[0].set_ylabel('$E_z$')
    axs[0].legend()
    axs[0].grid(True)

    # segundo subplot: H
    axs[1].plot(time_vec, Eplot,    label='Solução numérica')
    #axs[1].plot(time_vec, EplotEx,  label='Solução analítica')
    axs[1].set_title(f'Electrical field in ({xMon}, {yMon})')
    axs[1].set_xlabel('t(s)')
    axs[1].set_ylabel('$E_z$')
    axs[1].legend()
    axs[1].grid(True)

    plt.figure()
    #plt.plot(time_vec, Eplot2, label='Before Metasurface')
    plt.plot(time_vec, Eplot, label='After Metasurface')
    plt.xlabel('t(s)')
    plt.ylabel('E_z')
    plt.legend()
    plt.grid(True)

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

    return Hx, Hy, Ez, time

def plot_solution(x, y, u, label):

    # achata tudo para interpolação
    x_flat = x.flatten('F')
    y_flat = y.flatten('F')
    u_flat = u.flatten('F')

    # malha regular para contourf
    num_x, num_y = 200, 200
    xi = np.linspace(-a, a, num_x)
    yi = np.linspace(-b, b, num_y)
    Xi, Yi = np.meshgrid(xi, yi)

    Zi = griddata(
        points=(x_flat, y_flat),
        values=u_flat,
        xi=(Xi, Yi),
        method='cubic'
    )

    # --- plotting ---
    fig, ax = plt.subplots(figsize=(6,5))
    cf = ax.contourf(Xi, Yi, Zi, levels=50, cmap='viridis')
    plt.colorbar(cf, label=r'$E_z$')

    # 1) linha tracejada em x = xMeta
    ax.plot([xMeta, xMeta], [-yS, yS],
            color='red', linestyle='--', linewidth=2,
            label='metasurface')

    # 2) retângulo da fronteira TFSF
    tfsf_rect = Rectangle(
        (-xS, -yS),      # canto inferior esquerdo
        x0, y0,      # largura e altura
        fill=False,           # só contorno
        edgecolor='black',
        linestyle='-',
        linewidth=2,
        label='TFSF border'
    )
    ax.add_patch(tfsf_rect)

    # 3) retângulo da fronteira PML
    pml_rect = Rectangle(
        (-xpml, -ypml),
        2*xpml, 2*ypml,
        fill=False,
        edgecolor='cyan',
        linestyle='-',
        linewidth=2,
        label='PML border'
    )
    ax.add_patch(pml_rect)

    # ajustes finais
    ax.set_xlabel('$x$ (m)')
    ax.set_ylabel('$y$ (m)')
    ax.set_title(label)
    ax.axis('equal')

    # legenda
    ax.legend(loc='upper right')

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
        method='linear'
    )

    #faz o contourf no axes fornecido
    cf = ax.contourf(Xi, Yi, Zi, levels=50, vmin=vmin, vmax = vmax)
    ax.set_title(title)
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_aspect('equal')
    ax.set_xlim(-x0/2, x0/2)
    ax.set_ylim(-y0/2, y0/2)
    #Adiciona colorbal individual
    plt.colorbar(cf, ax=ax, label=r'$E_z$')

def ApplyPeriodicBDC2D(x, y, a, b, vM, vP, mM, mP):
    """
    Ajusta vmapP/mapP para impor CC periódica em x=±a e y=±b.

    Parâmetros:
      x, y      : arrays 1D de tamanho Np*K com coord. nodais (flattened columnwise)
      a, b      : meio-eixos do domínio em x e y
      NODETOL   : tolerância relativa para matching
      vmapM,P   : shape (Nfp, Nfaces, K) — antes de achatar
      mapM,P    : shape (Nfp, Nfaces, K) — indices lineares base para mapa
    Retorna:
      vmapP', mapP' ajustados
    """

    # coords de cada face-node
    xv = x.flatten('F')[vM]
    yv = y.flatten('F')[vM]

    # tolerâncias absolutas
    tolx = NODETOL * 2 * a
    toly = NODETOL * 2 * b

    # índices de fronteiras
    left   = np.where(np.abs(xv + a) < tolx)[0]
    right  = np.where(np.abs(xv - a) < tolx)[0]
    bottom = np.where(np.abs(yv + b) < toly)[0]
    top    = np.where(np.abs(yv - b) < toly)[0]

    # detectar cantos
    c_pp = np.where((np.abs(xv - a)   < tolx) & (np.abs(yv - b)   < toly))[0]  # (+a,+b)
    c_pm = np.where((np.abs(xv - a)   < tolx) & (np.abs(yv + b)   < toly))[0]  # (+a,-b)
    c_mp = np.where((np.abs(xv + a)   < tolx) & (np.abs(yv - b)   < toly))[0]  # (-a,+b)
    c_mm = np.where((np.abs(xv + a)   < tolx) & (np.abs(yv + b)   < toly))[0]  # (-a,-b)

    # remover cantos dos conjuntos unidimensionais
    def without(indices, remove):
        return np.setdiff1d(indices, remove, assume_unique=True)

    left   = without(left,   np.hstack([c_pp, c_pm, c_mp, c_mm]))
    right  = without(right,  np.hstack([c_pp, c_pm, c_mp, c_mm]))
    bottom = without(bottom, np.hstack([c_pp, c_pm, c_mp, c_mm]))
    top    = without(top,    np.hstack([c_pp, c_pm, c_mp, c_mm]))

    # função de casamento 1D
    def match(idA, idB, coordA, coordB, tol):
        for i in idA:
            diffs = np.abs(coordB[idB] - coordA[i])
            sel = np.where(diffs < tol)[0]
            if sel.size:
                j = idB[sel[0]]
                vP[i] = vM[j]
                mP[i] = mM[j]

    # casamentos em x e em y (sem cantos)
    match(left,  right,  yv, yv, toly)
    match(right, left,   yv, yv, toly)
    match(bottom, top,   xv, xv, tolx)
    match(top,    bottom,xv, xv, tolx)

    # agora tratar os 4 cantos: cada um anda diagonalmente 2a em x e 2b em y
    # helper para casar canto i -> canto oposto j
    def match_corner(ci, cj):
        if ci.size and cj.size:
            # assumimos correspondência 1 a 1
            for idx_i, idx_j in zip(ci, cj):
                vP[idx_i] = vM[idx_j]
                mP[idx_i] = mM[idx_j]

    # (+a,+b) ↔ (−a,−b)
    match_corner(c_pp, c_mm)
    match_corner(c_mm, c_pp)
    # (+a,−b) ↔ (−a,+b)
    match_corner(c_pm, c_mp)
    match_corner(c_mp, c_pm)

    return vP, mP

def maskMeta():
    #Função para calcular a máscara para a metasuperfície:

    # 1) decompõe vmapM em (linha, coluna)
    rows = vmapM % Np
    cols = vmapM // Np
    colsP =vmapP // Np

    # 2) pega coordenadas
    x_vals = x.flatten('F')[vmapM]
    y_vals = y.flatten('F')[vmapM]

    # 3) máscara conjunta
    maskX = (
        (np.isclose(x_vals, xMeta, rtol = NODETOL))
        & (np.abs(y_vals) <= yS)
    )

    # 4) contagem por coluna
    counts_per_colX = np.bincount(cols[maskX], minlength=K)

    # 5) colunas “válidas”
    valid_colsX = np.where(counts_per_colX > Nfp)[0]

    print("Colunas que têm ao menos", Nfp, 
        "pontos em x=xMeta:", valid_colsX)

    # 6) média de x em cada coluna (sobre todos os pontos do elemento)
    x_col_mean = x.mean(axis=0)   # vetor de tamanho K

    # 7) filtra colunas válidas de acordo com a média
    pre_cols_elemsX = valid_colsX[x_col_mean[valid_colsX] < xMeta]
    pos_cols_elemsX = valid_colsX[x_col_mean[valid_colsX] > xMeta]

    print("Colunas à esquerda (média x < -xS):",  pre_cols_elemsX)
    print("Colunas à direita (média x >  xS):",  pos_cols_elemsX)

    maskMetaPre = np.where(maskX & np.isin(cols, pre_cols_elemsX) & np.isin(colsP, pos_cols_elemsX))
    maskMetaPos = np.where(maskX & np.isin(cols, pos_cols_elemsX) & np.isin(colsP, pre_cols_elemsX))

    return maskMetaPre, maskMetaPos

def find_maskSource():
    # 1) decompõe vmapM em (linha, coluna)
    rows = vmapM % Np
    cols = vmapM // Np
    colsP =vmapP // Np

    # 2) pega coordenadas
    x_vals = x.flatten('F')[vmapM]
    y_vals = y.flatten('F')[vmapM]

    # 3) máscara conjunta
    #Se trocar o sinal de um dos xs impõe fonte sobre todo o domínio
    maskX = (
        (np.isclose(x_vals, -xS) | np.isclose(x_vals, xS))
        & (np.abs(y_vals) <= yS)
    )
    

    # 4) contagem por coluna
    counts_per_colX = np.bincount(cols[maskX], minlength=K)

    # 5) colunas “válidas”
    valid_colsX = np.where(counts_per_colX > Nfp)[0]

    print("Colunas que têm ao menos", Nfp, 
        "pontos em x=xS com |y|≤yS:", valid_colsX)

    # 6) média de x em cada coluna (sobre todos os pontos do elemento)
    x_col_mean = x.mean(axis=0)   # vetor de tamanho K

    # 7) filtra colunas válidas de acordo com a média
    sct_cols_elemsX  = valid_colsX[np.abs(x_col_mean[valid_colsX]) > xS]
    total_cols_elemsX = valid_colsX[np.abs(x_col_mean[valid_colsX]) < xS]

    print("Colunas em SFZ:",  sct_cols_elemsX)
    print("Colunas em TFS:",  total_cols_elemsX)

    #Para y
    # 3) máscara conjunta
    maskY = (
        (np.isclose(y_vals, yS) | np.isclose(y_vals, -yS))
        & (np.abs(x_vals) <= xS)
    )
    '''maskY = (
        (np.isclose(y_vals, yS) | np.isclose(y_vals, -yS))
        & (x_vals <= xMeta)
    )'''

    # 4) contagem por coluna
    counts_per_colY = np.bincount(cols[maskY], minlength=K)

    # 5) colunas “válidas”
    valid_colsY = np.where(counts_per_colY > Nfp)[0]

    print("Colunas que têm ao menos", Nfp, 
        "pontos em y=yS com |x|≤xS:", valid_colsY)

    # 6) média de x em cada coluna (sobre todos os pontos do elemento)
    y_col_mean = y.mean(axis=0)   # vetor de tamanho K

    # 7) filtra colunas válidas de acordo com a média
    sct_cols_elemsY  = valid_colsY[np.abs(y_col_mean[valid_colsY]) > yS]
    total_cols_elemsY = valid_colsY[np.abs(y_col_mean[valid_colsY]) < yS]

    print("Colunas em SFZ:",  sct_cols_elemsY)
    print("Colunas em TFZ:",  total_cols_elemsY)


    maskScatterRegion = np.where((maskX & np.isin(cols, sct_cols_elemsX) & np.isin(colsP, total_cols_elemsX)
                                  | (maskY & np.isin(cols, sct_cols_elemsY) & np.isin(colsP, total_cols_elemsY))))
    maskTotalRegion = np.where((maskX & np.isin(cols, total_cols_elemsX) & np.isin(colsP, sct_cols_elemsX)
                                  | (maskY & np.isin(cols, total_cols_elemsY) & np.isin(colsP, sct_cols_elemsY))))

    return maskScatterRegion, maskTotalRegion

def calculate_pml():
    #Definindo os parâmetros da PML
    #Condutância inicial
    sigma0x = 13.82
    sigma0y = 13.82
    p = 2

    #Calcula as coordenadas de x e y do centroide de cada elemento (vetores de tamanho K)
    x_col_mean = x.mean(axis=0)
    y_col_mean = y.mean(axis=0)

    #Inicializando os vetores
    sigma_x = np.zeros_like(x_col_mean)
    sigma_y = np.zeros_like(y_col_mean)

    #Em x
    maskXP = np.where(x_col_mean > xpml)
    maskXN = np.where(x_col_mean < -xpml)
    sigma_x[maskXP] = sigma0x * ((x_col_mean[maskXP] - xpml)**p)
    sigma_x[maskXN] = sigma0x * ((x_col_mean[maskXN] + xpml)**p)
    sigma_x = np.tile(sigma_x, (Np, 1)) #Replicar nos elementos

    #Em y
    maskYP = np.where(y_col_mean > ypml)
    maskYN = np.where(y_col_mean < -ypml)
    sigma_y[maskYP] = sigma0y * ((y_col_mean[maskYP] - ypml)**p)
    sigma_y[maskYN] = sigma0y * ((y_col_mean[maskYN] + ypml)**p)
    sigma_y = np.tile(sigma_y, (Np, 1)) #Replicar nos elementos
    
    return sigma_x, sigma_y

#Definições iniciais
N = 3
Np = int((N+1) * (N+2) / 2)
Nfp = N +1
Nfaces = 3
NODETOL = 1E-15
eps = 1e-16
alpha = 1

#Variáveis de controle
metaOn = 0 #metaOn = 1 -> metasuperfície ligada
SourceTFSF = 1 
SourceGaussian = 0

#Pontos no triângulo equilátero
xeq, yeq = Nodes2D(N)

#Pontos no triângulo de referência
r, s = xytors(xeq, yeq)

#Matrizes necessárias
V = Vandermonde2D(N, r, s)
Dr, Ds = Dmatrices2D(N, r, s, V)
invV = np.linalg.inv(V)
MassMatrix = invV.T @ invV

#Criando a malha
#Parâmetros físicos
L = 1 #Comprimento da PML
x0 = 2 #Tamanho do domínio físico em x (-x0/2 a x0/2)
y0 = 2 #Tamanho do domínio físico em y (-y0/2 a y0/2)
LS = 0 #Comprimento do espaço TFSF
xpml = x0/2 + LS #Coordenada x de início da PML
ypml = y0/2 + LS #Coordenada y de início da PML
xS = x0/2
yS = y0/2
a = L + LS + x0/2 #Comprimento máximo do eixo x a partir da origem
b = L + LS + y0/2 #Comprimento máximo do eixo y a partir da origem

#Definindo a metasuperfície, paralela ao plano y-z
xMeta = 0
T = 1
R = np.sqrt(1 - T**2)
#R = 0

#Gerando a malha
num_elements = 16
vertices, EToV = generate_mesh(num_elements, a, b, x0, y0, xpml, ypml, xMeta)
VX = vertices[:, 0]
VY = vertices[:, 1]

VX = np.load("C:\\git\\PyDG1D\\examplesData\\outputs\\cem_6\\SquareK16_pyData\\vx.npy")
VY = np.load("C:\\git\\PyDG1D\\examplesData\\outputs\\cem_6\\SquareK16_pyData\\vy.npy")
vertices = np.column_stack((VX, VY))
EToV = np.load("C:\\git\\PyDG1D\\examplesData\\outputs\\cem_6\\SquareK16_pyData\\EToV.npy")

plot_mesh(vertices, EToV, order=1)

Nv = len(VX)
K = len(EToV) #Número de elementos
print(f'número de elementos: {K}')

#Encontrando os pontos dentro do triângulo:
va = EToV[:, 0]
vb = EToV[:, 1]
vc = EToV[:, 2]
x = 0.5*(-(r+s)[:, None] @ VX[None, va] + (1+r)[:, None] @ VX[None, vb] + (1+s)[:, None] @ VX[None, vc])
y = 0.5*(-(r+s)[:, None] @ VY[None, va] + (1+r)[:, None] @ VY[None, vb] + (1+s)[:, None] @ VY[None, vc])

#Encontrando os pontos que estão nas interfaces (arestas) dos triângulos
fmask1 = np.where(np.abs(s + 1) < NODETOL)
fmask2 = np.where(np.abs(r + s) < NODETOL)
fmask3 = np.where(np.abs(r + 1) < NODETOL)
Fmask = np.concatenate([fmask1, fmask2, fmask3]).T
Fx = x[Fmask.flatten('F')]
Fy = y[Fmask.flatten('F')]

#Matriz de Lifting
LIFT = Lift2D()

#Fatores geométricos
rx, sx, ry, sy, J = GeometricFactors2D(x, y, Dr, Ds)
nx, ny, sJ = Normals2D() #vetores normais
Fscale = sJ/J[Fmask.flatten('F')]


#Matrizes e vetores de mapeamento
EToE, EToF = Connect2D(EToV)
mapM, mapP, vmapM, vmapP, vmapB, mapB = BuildMaps2D()

#Calcula os operadores fracos
Vr, Vs = GradVandermonde2D(N, r, s)
Drw = np.linalg.solve((V @ V.T), (V @ Vr.T).T).T
Dsw = np.linalg.solve((V @ V.T), (V @ Vs.T).T).T

#Coeficientes de runge-kutta
rk4a, rk4b, rk4c = rungekutta()

sigma_x, sigma_y = calculate_pml()

if metaOn == 1:
    #Calculando a máscara para metasuperfície
    maskMetaleft, maskMetaright = maskMeta()
    print(f'Máscara meta a esquerda:{maskMetaleft}')
    print(f'Máscara meta a direita:{maskMetaright}')


#Calculando máscara para fonte
if SourceTFSF == 1:
    maskScatterRegion, maskTotalRegion = find_maskSource()
    '''print(f'Máscara da fonte a esquerda: {len(x.flatten('F')[vmapM[maskScatterRegion]])}')
    print(f'Máscara fonte da direita: {len(x.flatten('F')[vmapM[maskTotalRegion]])}')'''

    #Inicializando o problema
    #Condição inicial
    Ez = np.zeros((Np, K))
    Hx = np.zeros((Np, K))
    Hy = np.zeros((Np, K))
    #Sinais auxiliares da PML
    Px = np.zeros((Np, K))
    Py = np.zeros((Np, K))
    Pz = np.zeros((Np, K))
    #Onda incidente
    #Comprimento de onda
    E0 = 1
    lambda_ = x0/4
    k = 2 * np.pi / lambda_
    f = 1/lambda_
    omega = 2 * np.pi * f

if SourceGaussian == 1:
    #Inicializando o problema
    xFonte = 0
    yFonte = 0
    stddv = (x0/2)/4
    #Condição inicial
    Ez = np.exp(-(x - xFonte)**2/(2*stddv**2)) * np.exp(-(y - yFonte)**2/(2*stddv**2))
    Hx = np.zeros((Np, K))
    Hy = np.zeros((Np, K))
    #Sinais auxiliares da PML
    Px = np.zeros((Np, K))
    Py = np.zeros((Np, K))
    Pz = np.zeros((Np, K))

plot_solution(x, y, Ez, 'Condição inicial')

#Define e aplica condições de contorno 
#Contour = 0 -> periódicas
#Contour = 1 -> PEC
#Outros valores -> Sem condições de contorno.
contour = 1
if contour == 0:
    vmapP, mapP = ApplyPeriodicBDC2D(x, y, a, b, vmapM, vmapP, mapM, mapP)

#Tempo final de simulação
FinalTime = 1

#Chama a rotina principal de integração temporal
Hx, Hy, Ez, time = Maxwell2D(Hx, Hy, Ez, Px, Py, Pz, FinalTime)

#Impressão de confirmação
print(f'Simulação concluída em t = {time:.3f}')

plot_solution(x, y, Ez, '')
#plot_solution(x, y, Ez_ex, 'Solução analítica')

DIR_PATH = "C:\\git\\PyDG1D\\examplesData\\outputs\\cem_6\\SquareK16_pyData\\"

vmapI_tfz_scatter = vmapM[maskScatterRegion]
vmapI_tfz_total = vmapM[maskTotalRegion]

# print(f'vmapM: {vmapM.shape}: \n{vmapM}')
# print(f'vmapP: {vmapP.shape}: \n{vmapP}')

assert np.array_equal(EToE, np.load(DIR_PATH + "EToE.npy")), "EToE"
assert np.array_equal(EToF, np.load(DIR_PATH + "EToF.npy")), "EToF"

assert np.array_equal(vmapM, np.load(DIR_PATH + "vmapM.npy")), "vmapM"
assert np.array_equal(vmapP, np.load(DIR_PATH + "vmapP.npy")), "vmapP"
assert np.array_equal(vmapB, np.load(DIR_PATH + "vmapB.npy")), "vmapB"
assert np.array_equal(mapB, np.load(DIR_PATH + "mapB.npy")), "mapB"

assert np.array_equal(vmapI_tfz_scatter, np.load(DIR_PATH + "vmapI_tfz_scatter.npy")), "vmapI_tfz_scatter"
#assert np.array_equal(vmapI_tfz_total, np.load(DIR_PATH + "vmapI_tfz_total.npy")), "vmapI_tfz_total"


# plt.show()

