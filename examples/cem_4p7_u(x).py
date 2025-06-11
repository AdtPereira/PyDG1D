import os
import sys

# Caminho do projeto
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.tri import Triangulation

from maxwell.dg.dg2d import *
from maxwell.utils import *
from maxwell.dg.dg2d_tools import *
from mesher.create_mesh import *


BOUNDARY = [{'tag': 101,
                 'type': 'Dirichlet',
                 'value': 0.0,
                 'name': 'Ez_0'}]
       
MATERIAL = [{'tag': 201,
                 'name': 'free_space',
                 'relative_magnetic_permeability': 1,
                 'relative_electric_permittivity': 1}]  


def exact_u(x, y):
    """Função exata para análise de convergência."""
    return np.cos(np.pi * x) * np.cos(np.pi * y)


def build_reference_interpolation(N, exact_func):
    """
    Retorna os coeficientes da interpolação da função exata no triângulo de referência.

    Retorna:
    --------
    coeff : ndarray
        Coeficientes de interpolação da função exata nos nós nodais.
    V : ndarray
        Matriz de Vandermonde nos nós nodais (para reuso).
    """
    x, y = set_nodes_in_equilateral_triangle(N)
    r, s = xy_to_rs(x, y)
    V = vandermonde(N, r, s)
    u_nodes = exact_func(r, s)
    coeff = np.linalg.inv(V) @ u_nodes
    return coeff, V


def compute_L2_error(u_exact, u_approx, wq):
    """Calcula a norma L2 discreta do erro."""
    diff = u_exact - u_approx
    return np.sqrt(np.sum(wq * diff**2))


def compute_error_L2_reference_I(N):
    """
    Calcula o erro L2 da interpolação da função exata no triângulo de referência.
    """
    coeff, V = build_reference_interpolation(N, exact_u)

    # Pontos de quadratura
    Nq = N + 3
    xq, yq = set_nodes_in_equilateral_triangle(Nq)
    rq, sq = xy_to_rs(xq, yq)
    wq = np.full_like(rq, 2.0 / len(rq))  # área total do triângulo de referência
    Vq = vandermonde(N, rq, sq)

    u_interp = Vq @ coeff
    u_exact = exact_u(rq, sq)

    return compute_L2_error(u_exact, u_interp, wq)


def compute_error_L2_global(N, mesh_data, exact_func):
    """
    Calcula o erro global L2 no domínio físico usando a interpolação polinomial por elemento.
    """
    VX, VY, EToV = mesh_data['VX'], mesh_data['VY'], mesh_data['EToV']
    K = EToV.shape[0]
    
    coeff, V = build_reference_interpolation(N, exact_func)

    # Quadratura no triângulo de referência
    Nq = N + 3
    xq, yq = set_nodes_in_equilateral_triangle(Nq)
    rq, sq = xy_to_rs(xq, yq)
    Vq = vandermonde(N, rq, sq)
    u_interp_local = Vq @ coeff
    wq = np.full_like(rq, 2.0 / len(rq))

    error_total = 0.0

    for k in range(K):
        va, vb, vc = EToV[k]
        xk = 0.5 * (-(rq + sq) * VX[va] + (1 + rq) * VX[vb] + (1 + sq) * VX[vc])
        yk = 0.5 * (-(rq + sq) * VY[va] + (1 + rq) * VY[vb] + (1 + sq) * VY[vc])

        Jk = abs((VX[vb] - VX[va]) * (VY[vc] - VY[va]) -
                 (VX[vc] - VX[va]) * (VY[vb] - VY[va])) / 2

        u_exact_k = exact_func(xk, yk)
        diff2 = (u_exact_k - u_interp_local) ** 2
        error_total += np.sum(wq * diff2) * Jk

    return np.sqrt(error_total)


def build_convergence_table(N_list, errors):
    """
    Constrói um DataFrame com N, erro L2 e taxa de convergência.

    Parâmetros
    ----------
    N_list : list[int]
        Lista das ordens dos polinômios.
    errors : list[float]
        Lista dos erros L2 associados a cada ordem.

    Retorno
    -------
    df : pd.DataFrame
        DataFrame contendo colunas: Ordem N, Erro L2, Taxa de Convergência
    """
    rates = [np.nan]  # sem taxa para o primeiro ponto
    for i in range(1, len(N_list)):
        e0, e1 = errors[i - 1], errors[i]
        N0, N1 = N_list[i - 1], N_list[i]
        if e0 > 0 and e1 > 0:
            rate = np.log(e0 / e1) / np.log(N1 / N0)
            rate = round(rate, 2)  # arredondar para 2 casas decimais
        else:
            rate = np.nan
        rates.append(rate)

    df = pd.DataFrame({
        "Ordem N": N_list,
        "Erro L2": errors,
        "Taxa de Convergência": rates
    })
    return df


def plot_convergence(N_list, errors, label="Erro $L^2$", title="", filename=None):
    plt.figure(figsize=(8, 5))
    plt.semilogy(N_list, errors, 'o-', label=label)

    # --- Comparação com curva analítica esperada ---
    alpha = 1.0
    C = errors[0] * np.exp(alpha * N_list[0])  # normalização
    expected = C * np.exp(-alpha * np.array(N_list))
    plt.semilogy(N_list, expected, '--', label=f"$Ce^{{-{alpha}N}}$")

    # --- Estética ---
    plt.xlabel("Ordem N")
    plt.ylabel("Erro")
    plt.title(title or "Convergência espectral")
    plt.grid(True, which='both', linestyle='--', linewidth=0.5)
    plt.legend()
    plt.tight_layout()

    if filename:
        plt.savefig(filename, dpi=300)
    else:
        plt.show()


def interpolate_solution_on_physical_domain(N, mesh_data, exact_func):
    """
    Interpola a função exata u(x, y) nos pontos de colocação do domínio físico,
    e retorna os valores para visualização.

    Parâmetros:
    -----------
    N : int
        Ordem do polinômio.
    mesh_data : dict
        Deve conter 'VX', 'VY' e 'EToV'.
    exact_func : callable
        Função exata u(x, y) a ser interpolada.

    Retorna:
    --------
    x_plot, y_plot : (Np_plot, K)
        Coordenadas físicas dos pontos de visualização.
    u_plot : (Np_plot,)
        Valores interpolados da função nos pontos (x_plot, y_plot).
    tri_ref : matplotlib.tri.Triangulation
        Malha de referência usada para renderização.
    """

    # 1. Pontos nodais e matriz de interpolação no triângulo de referência
    coeff, V = build_reference_interpolation(N, exact_func)

    # 2. Pontos de visualização mais densos no triângulo de referência
    N_plot = max(N + 3, 6)
    r_plot, s_plot = xy_to_rs(*set_nodes_in_equilateral_triangle(N_plot))
    V_plot = vandermonde(N, r_plot, s_plot)
    u_plot = V_plot @ coeff

    # 3. Coordenadas físicas dos pontos (r_plot, s_plot)
    x_nodes, y_nodes = nodes_coordinates_from_gmsh(N_plot, mesh_data)

    # 4. Triangulação de referência (a mesma para todos os elementos)
    tri_ref = Triangulation(r_plot, s_plot)

    interpolate_data = {
        'x_plot': x_nodes,
        'y_plot': y_nodes,
        'u_plot': u_plot,   
        'tri_ref': tri_ref
    }

    return interpolate_data


def plot_collocation_points(N, mesh_data, show_ids=False):
    """
    Plota os pontos de colocação no domínio físico, com marcador distinto por elemento.
    Os índices são desenhados deslocados em direção ao centróide do elemento.

    Parâmetros:
    -----------
    N : int
        Ordem polinomial dos elementos.
    mesh_data : dict
        Dicionário com 'VX', 'VY' e 'EToV'.
    show_ids : bool
        Se True, mostra os índices dos nós locais deslocados para evitar sobreposição.
    """

    VX = mesh_data['VX']
    VY = mesh_data['VY']
    EToV = mesh_data['EToV']
    K = EToV.shape[0]

    x_nodes, y_nodes = nodes_coordinates_from_gmsh(N, mesh_data)
    markers = ['o', 's', '^', 'D', 'v', '>', '<', 'p', '*', 'h', 'X']
    colors = plt.cm.tab10(np.linspace(0, 1, 10))

    # Malha de triângulos
    for k in range(K):
        nodes = EToV[k]
        poly_x = VX[nodes]
        poly_y = VY[nodes]
        plt.plot(np.append(poly_x, poly_x[0]), np.append(poly_y, poly_y[0]), 'k-', lw=0.8)

    # Pontos de colocação com marcador distinto por elemento
    for k in range(K):
        xk = x_nodes[:, k]
        yk = y_nodes[:, k]
        mk = markers[k % len(markers)]
        ck = colors[k % len(colors)]
        plt.scatter(xk, yk, marker=mk, color=ck, s=35, label=f'e{k}')

        if show_ids:
            # centro do elemento
            va, vb, vc = EToV[k]
            xc = (VX[va] + VX[vb] + VX[vc]) / 3
            yc = (VY[va] + VY[vb] + VY[vc]) / 3

            for i in range(x_nodes.shape[0]):
                # vetor do ponto para o centro
                dx = xc - xk[i]
                dy = yc - yk[i]
                scale = 0.1  # fator de deslocamento
                xt = xk[i] + dx * scale
                yt = yk[i] + dy * scale
                plt.text(xt, yt, str(i), fontsize=7, color='black', ha='center', va='center')

    plt.gca().set_aspect('equal')
    plt.xlabel("x")
    plt.ylabel("y")
    plt.title(fr"Pontos de colocação no domínio físico ($N={N}$)")
    plt.grid(True, linestyle='--', linewidth=0.3)
    plt.legend(fontsize=8, loc='center left', bbox_to_anchor=(1, 0.5))
    plt.tight_layout()
    plt.show()


def plot_collocation_and_interpolated_points(interpolate_data, N, mesh_data):
    """
    Compara visualmente os pontos de colocação com os pontos de visualização interpolados.

    Parâmetros:
    -----------
    N : int
        Ordem do polinômio (usada para os pontos de colocação).
    mesh_data : dict
        Dicionário contendo 'VX', 'VY', 'EToV'.
    x_plot, y_plot : ndarray (Np_plot, K)
        Pontos físicos de visualização interpolados, retornados por interpolate_solution_on_physical_domain.
    """
    # Pontos de colocação nos elementos
    x_plot, y_plot = interpolate_data['x_plot'], interpolate_data['y_plot']
    x_nodes, y_nodes = nodes_coordinates_from_gmsh(N, mesh_data)
    K = x_nodes.shape[1]

    # Plotagem
    plt.figure(figsize=(8, 6))

    # 1. Pontos de colocação
    for k in range(K):
        plt.scatter(x_nodes[:, k], y_nodes[:, k],
                    marker='o', color='blue', s=35, label='Colocação' if k == 0 else "")

    # 2. Pontos interpolados (visualização)
    for k in range(K):
        plt.scatter(x_plot[:, k], y_plot[:, k],
                    marker='x', color='red', s=25, label='Interpolação' if k == 0 else "")

    # 3. Opcional: bordas dos elementos
    VX = mesh_data['VX']
    VY = mesh_data['VY']
    EToV = mesh_data['EToV']
    for k in range(K):
        nodes = EToV[k]
        poly_x = VX[nodes]
        poly_y = VY[nodes]
        plt.plot(np.append(poly_x, poly_x[0]), np.append(poly_y, poly_y[0]), 'k-', lw=0.8)

    plt.xlabel("x")
    plt.ylabel("y")
    plt.title(f"Comparação entre pontos de colocação e interpolação (N = {N})")
    plt.gca().set_aspect('equal')
    plt.grid(True, linestyle='--', linewidth=0.4)
    plt.legend()
    plt.tight_layout()
    plt.show()


def analyze_convergence(N_list, mesh_data, exact_func, compute_error_fn, label, title):
    """
    Realiza a análise de convergência para uma lista de ordens N.

    Parâmetros:
    -----------
    N_list : iterable[int]
        Lista de ordens polinomiais.
    mesh_data : dict or None
        Malha usada no cálculo. Pode ser None para erro no triângulo de referência.
    exact_func : callable
        Função exata u(x, y).
    compute_error_fn : callable
        Função que recebe (N, ...) e retorna o erro L2.
    label : str
        Rótulo da curva no gráfico.
    title : str
        Título do gráfico.
    """
    errors = [compute_error_fn(N, mesh_data, exact_func) if mesh_data else compute_error_fn(N)
              for N in N_list]
    
    df = build_convergence_table(list(N_list), errors)
    print(f"\n{title}")
    print(df.to_string(index=False))

    plot_convergence(N_list, errors, label=label, title=title)


def build_global_visualization_mesh(N_plot, mesh_data):
    """
    Gera uma malha global de visualização a partir da malha física, com conectividades contínuas.

    Retorna:
    --------
    x_all, y_all : ndarray
        Todos os pontos físicos (flattened) da visualização.
    triangles_all : ndarray
        Conectividade global dos triângulos (para Triangulation).
    """
    VX, VY, EToV = mesh_data['VX'], mesh_data['VY'], mesh_data['EToV']
    K = EToV.shape[0]

    x_ref, y_ref = set_nodes_in_equilateral_triangle(N_plot)
    r_plot, s_plot = xy_to_rs(x_ref, y_ref)
    tri_ref = Triangulation(r_plot, s_plot)
    tri_local = tri_ref.triangles
    Np = len(r_plot)

    x_all = []
    y_all = []
    triangles_all = []
    node_offset = 0

    for k in range(K):
        va, vb, vc = EToV[k]
        xk = 0.5 * (-(r_plot + s_plot) * VX[va] + (1 + r_plot) * VX[vb] + (1 + s_plot) * VX[vc])
        yk = 0.5 * (-(r_plot + s_plot) * VY[va] + (1 + r_plot) * VY[vb] + (1 + s_plot) * VY[vc])

        x_all.extend(xk)
        y_all.extend(yk)
        triangles_all.extend(tri_local + node_offset)
        node_offset += Np

    return np.array(x_all), np.array(y_all), np.array(triangles_all)


def plot_exact_solution_global(N_plot, mesh_data, exact_func, field_label='u'):
    """
    Plota a solução exata u(x,y) sobre uma malha global contínua.
    """
    x_all, y_all, triangles_all = build_global_visualization_mesh(N_plot, mesh_data)
    u_exact = exact_func(x_all, y_all)

    # --- Gráfico 2D e 3D lado a lado ---
    fig = plt.figure(figsize=(14, 6))
    ax1 = fig.add_subplot(1, 2, 1)
    ax2 = fig.add_subplot(1, 2, 2, projection='3d')

    tri_global = Triangulation(x_all, y_all, triangles_all)

    # Contorno 2D
    tcf = ax1.tricontourf(tri_global, u_exact, cmap='viridis', levels=100)
    fig.colorbar(tcf, ax=ax1, label=fr'{field_label}(x,y) exato')
    ax1.set_title(fr'Solução exata ${field_label}(x, y)$ - Contorno Global')
    ax1.set_xlabel("x")
    ax1.set_ylabel("y")
    ax1.set_aspect("equal")
    ax1.grid(True, linestyle='--', linewidth=0.3)

    # Superfície 3D
    ax2.plot_trisurf(tri_global, u_exact, cmap='viridis', linewidth=0.2)
    ax2.set_title(fr'Solução exata ${field_label}(x, y)$ - Superfície')
    ax2.set_xlabel("x")
    ax2.set_ylabel("y")
    ax2.set_zlabel(fr'{field_label}(x, y)')
    ax2.view_init(elev=35, azim=-120)

    plt.tight_layout()
    plt.show()


def main():
    clear_terminal()
    PROBLEM = {'name': 'cem_4p7_u(x)',
        'folder_name': 'cem_4p7_u(x)',
        'description': '',
        'L': 2.0,
        'n_order': 3
    }

    # Criar a malha retangular com Gmsh
    mesh_data = mesh_rectangular_domain(PROBLEM, BOUNDARY, MATERIAL, h=2, view_mesh=False, mesh_info=False)

    interpolate_data = interpolate_solution_on_physical_domain(
        PROBLEM['n_order'], mesh_data, exact_func=exact_u)
    
    # Plotar os pontos de colocação (nós nodais) no domínio físico
    plot_collocation_points(PROBLEM['n_order'], mesh_data, show_ids=False)

    # Comparar pontos de colocação com pontos interpolados
    plot_collocation_and_interpolated_points(interpolate_data, PROBLEM['n_order'], mesh_data)

    # Plotar solução exata no domínio físico
    plot_exact_solution_global(
        N_plot=PROBLEM['n_order'] + 3,
        mesh_data=mesh_data,
        exact_func=exact_u,
        field_label='u'
    )

    # Erro no triângulo de referência
    N_list = range(1, 15)
    analyze_convergence(
        N_list=N_list,
        mesh_data=None,
        exact_func=exact_u,
        compute_error_fn=compute_error_L2_reference_I,
        label="Erro $L^2$",
        title="Convergência espectral no triângulo de referência"
    )

    


if __name__ == '__main__':
    main()
