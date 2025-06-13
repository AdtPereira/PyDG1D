import os
import sys

# Caminho do projeto
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..\..')))

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.tri import Triangulation

from maxwell.dg.dg2d import *
from maxwell.utils import *
from maxwell.dg.dg2d_tools import *
from mesher.create_mesh import *

BOUNDARY = [{'tag': 101, 'type': 'Dirichlet', 'value': 0.0, 'name': 'Ez_0'}]
MATERIAL = [{'tag': 201, 'name': 'free_space', 'relative_magnetic_permeability': 1, 'relative_electric_permittivity': 1}]
PROBLEM = {'name': 'cem_4p7', 'folder_name': 'cem_4p7', 'L': 2.0}


class SpectralAnalyzer:
    def __init__(self, mesh, n_order, exact_func):
        self.mesh = mesh
        self.exact_func = exact_func
        self.N = n_order                    # Ordem polinomial dos elementos
        self.N_plot = self.N + 2            # Ordem dos pontos de visualização


    def build_reference_interpolation(self):
        """
        Gera os nós de colocação no triângulo de referência e constrói a matriz de Vandermonde  
        para interpolação de uma função exata definida no triângulo de referência.
        """
        x, y = set_nodes_in_equilateral_triangle(self.N)
        r, s = xy_to_rs(x, y)
        V = vandermonde(self.N, r, s)
        u_nodes = self.exact_func(r, s)
        coeff = np.linalg.inv(V) @ u_nodes
        return coeff, V
    

    def interpolate_solution(self):
        """
        Interpola a solução exata u(x,y) para pontos físicos de visualização
        distribuídos em todos os elementos da malha.

        Retorna:
            dict contendo 'x_plot', 'y_plot', 'u_plot' com shape (Np_plot, K)
        """
        EToV = self.mesh.EToV
        K = EToV.shape[0]

        # Pontos de visualização no triângulo de referência
        x_ref, y_ref = set_nodes_in_equilateral_triangle(self.N + 1)
        r_plot, s_plot = xy_to_rs(x_ref, y_ref)
        V_plot = vandermonde(self.N, r_plot, s_plot)
        Np_plot = len(r_plot)

        # Coordenadas físicas dos pontos de visualização e nodais
        x_plot, y_plot = nodes_coordinates(self.N + 1, self.mesh)
        x_nodes, y_nodes = nodes_coordinates(self.N, self.mesh)

        u_plot = np.zeros((Np_plot, K))
        r_nodes, s_nodes = xy_to_rs(*set_nodes_in_equilateral_triangle(self.N))
        V = vandermonde(self.N, r_nodes, s_nodes)

        for k in range(K):
            u_nodes = self.exact_func(x_nodes[:, k], y_nodes[:, k])
            coeff_k = np.linalg.solve(V, u_nodes)
            u_plot[:, k] = V_plot @ coeff_k

        return {
            'x_plot': x_plot,
            'y_plot': y_plot,
            'u_plot': u_plot,
            'tri_ref': Triangulation(r_plot, s_plot)
        }


    def L2_error(self, u_exact, u_approx, wq):
        return np.sqrt(np.sum(wq * (u_exact - u_approx)**2))


    def reference_element_error(self):
        coeff, V = self.build_reference_interpolation()
        xq, yq = set_nodes_in_equilateral_triangle(self.N_plot)
        rq, sq = xy_to_rs(xq, yq)
        Vq = vandermonde(self.N, rq, sq)
        u_interp = Vq @ coeff
        u_exact = self.exact_func(rq, sq)
        wq = np.full_like(rq, 2.0 / len(rq))
        return self.L2_error(u_exact, u_interp, wq)


    def convergence_table(self, N_list, errors):
        rates = [np.nan]
        for i in range(1, len(N_list)):
            e0, e1 = errors[i - 1], errors[i]
            rate = np.log(e0 / e1) / np.log(N_list[i] / N_list[i - 1]) if e0 > 0 and e1 > 0 else np.nan
            rates.append(round(rate, 2))
        return pd.DataFrame({"N": N_list, "Erro L2": errors, "Taxa": rates})


    def global_domain_error(self):
        """
        Calcula o erro L2 global da solução interpolada no domínio físico.

        Retorna:
            erro_L2 : float
                Norma L2 do erro entre solução exata e interpolada no domínio físico.
        """
        data = self.interpolate_solution()

        x_plot = data['x_plot']
        y_plot = data['y_plot']
        u_interp = data['u_plot']
        u_exact = self.exact_func(x_plot, y_plot)

        Np_plot = x_plot.shape[0]
        wq = np.full(Np_plot, 2.0 / Np_plot)  # peso médio por ponto (área total ≈ 2)

        # Soma dos erros quadráticos por elemento
        err_quad = wq @ ((u_exact - u_interp) ** 2)
        return np.sqrt(np.sum(err_quad))


class SpectralVisualizer():
    def __init__(self, analyzer, mesh):
        self.mesh = mesh
        self.N = analyzer.N             # Ordem polinomial dos elementos
        self.N_plot = analyzer.N_plot   # Ordem dos pontos de visualização


    def build_global_visualization_mesh(self):
        """
        Gera uma malha global de visualização a partir da malha física, com conectividades contínuas.

        Retorna:
        --------
        x_all, y_all : ndarray
            Todos os pontos físicos (flattened) da visualização.
        triangles_all : ndarray
            Conectividade global dos triângulos (para Triangulation).
        """
        EToV = self.mesh.EToV
        K = EToV.shape[0]

        # Coordenadas no triângulo de referência
        x_ref, y_ref = set_nodes_in_equilateral_triangle(self.N_plot)
        r_plot, s_plot = xy_to_rs(x_ref, y_ref)
        tri_ref = Triangulation(r_plot, s_plot)
        tri_local = tri_ref.triangles
        Np = len(r_plot)

        # Coordenadas físicas globais dos pontos de visualização
        x_plot, y_plot = nodes_coordinates(self.N_plot, self.mesh)

        x_all = []
        y_all = []
        triangles_all = []
        node_offset = 0

        for k in range(K):
            x_all.extend(x_plot[:, k])
            y_all.extend(y_plot[:, k])
            triangles_all.extend(tri_local + node_offset)
            node_offset += Np

        return np.array(x_all), np.array(y_all), np.array(triangles_all)


    def collocation_points(self, show_ids=False):
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
        VX, VY, EToV = self.mesh.vx, self.mesh.vy, self.mesh.EToV
        K = EToV.shape[0]

        x_nodes, y_nodes = nodes_coordinates(self.N, self.mesh)
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
            plt.scatter(xk, yk, marker=mk, facecolors='none', color=ck, s=35, label=f'e{k}')

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
        plt.title(fr"Pontos de colocação no domínio físico ($N={self.N}$)")
        plt.grid(True, linestyle='--', linewidth=0.3)
        plt.legend(fontsize=10, loc='center left', bbox_to_anchor=(1, 0.5))
        plt.tight_layout()
        plt.show()


    def collocation_and_interpolated_points(self, interpolate_data):
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
        x_nodes, y_nodes = nodes_coordinates(self.N, self.mesh)
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
        VX, VY, EToV = self.mesh.vx, self.mesh.vy, self.mesh.EToV
        for k in range(K):
            nodes = EToV[k]
            poly_x = VX[nodes]
            poly_y = VY[nodes]
            plt.plot(np.append(poly_x, poly_x[0]), np.append(poly_y, poly_y[0]), 'k-', lw=0.8)

        plt.xlabel("x")
        plt.ylabel("y")
        plt.title(f"Comparação entre pontos de colocação e interpolação (N = {self.N})")
        plt.gca().set_aspect('equal')
        plt.grid(True, linestyle='--', linewidth=0.4)
        plt.legend(fontsize=10, loc='center left', bbox_to_anchor=(1, 0.5))
        plt.tight_layout()
        plt.show()


    def exact_solution_global(self, exact_func, field_label='u'):
        """
        Plota a solução exata u(x,y) sobre uma malha global contínua.
        """
        x_all, y_all, triangles_all = self.build_global_visualization_mesh()
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


    def interpolated_global_solution(self, analyzer, field_label='u'):
        """
        Plota a solução interpolada u(x,y) globalmente com continuidade entre elementos.
        Utiliza avaliação local por elemento e unifica a malha e a solução interpolada.
        """
        N = self.N
        VX, VY, EToV = self.mesh.vx, self.mesh.vy, self.mesh.EToV
        K = EToV.shape[0]

        # Pontos de visualização no triângulo de referência
        x_ref, y_ref = set_nodes_in_equilateral_triangle(N + 1)
        r_plot, s_plot = xy_to_rs(x_ref, y_ref)
        V_plot = vandermonde(N, r_plot, s_plot)
        tri_ref = Triangulation(r_plot, s_plot)
        tri_local = tri_ref.triangles
        Np_plot = len(r_plot)

        x_all, y_all, u_all = [], [], []
        triangles_all = []
        node_offset = 0

        for k in range(K):
            va, vb, vc = EToV[k]
            xk = 0.5 * (-(r_plot + s_plot) * VX[va] + (1 + r_plot) * VX[vb] + (1 + s_plot) * VX[vc])
            yk = 0.5 * (-(r_plot + s_plot) * VY[va] + (1 + r_plot) * VY[vb] + (1 + s_plot) * VY[vc])

            # Avaliação da função exata nos nós do elemento
            r_nodes, s_nodes = xy_to_rs(*set_nodes_in_equilateral_triangle(N))
            x_nodes = 0.5 * (-(r_nodes + s_nodes) * VX[va] + (1 + r_nodes) * VX[vb] + (1 + s_nodes) * VX[vc])
            y_nodes = 0.5 * (-(r_nodes + s_nodes) * VY[va] + (1 + r_nodes) * VY[vb] + (1 + s_nodes) * VY[vc])
            u_nodes = analyzer.exact_func(x_nodes, y_nodes)

            V = vandermonde(N, r_nodes, s_nodes)
            coeff_k = np.linalg.solve(V, u_nodes)
            u_k = V_plot @ coeff_k

            x_all.extend(xk)
            y_all.extend(yk)
            u_all.extend(u_k)
            triangles_all.extend(tri_local + node_offset)
            node_offset += Np_plot

        x_all = np.array(x_all)
        y_all = np.array(y_all)
        u_all = np.array(u_all)
        triangles_all = np.array(triangles_all)
        tri_global = Triangulation(x_all, y_all, triangles_all)

        # Gráficos 2D e 3D
        fig = plt.figure(figsize=(14, 6))
        ax1 = fig.add_subplot(1, 2, 1)
        ax2 = fig.add_subplot(1, 2, 2, projection='3d')

        tcf = ax1.tricontourf(tri_global, u_all, levels=100, cmap='viridis')
        fig.colorbar(tcf, ax=ax1, label=fr'{field_label}(x,y) interpolado')
        ax1.set_title(fr'Solução interpolada ${field_label}(x, y)$ - Contorno Global')
        ax1.set_xlabel("x")
        ax1.set_ylabel("y")
        ax1.set_aspect("equal")
        ax1.grid(True, linestyle='--', linewidth=0.3)

        ax2.plot_trisurf(tri_global, u_all, cmap='viridis', linewidth=0.2)
        ax2.set_title(fr'Solução interpolada ${field_label}(x, y)$ - Superfície')
        ax2.set_xlabel("x")
        ax2.set_ylabel("y")
        ax2.set_zlabel(fr'{field_label}(x, y)')
        ax2.view_init(elev=35, azim=-120)

        plt.tight_layout()
        plt.show()


    def convergence_rate(self, N_list, errors, label="Erro $L^2$", title="", filename=None):
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


def exact_u(x, y):
    return np.cos(np.pi * x) * np.cos(np.pi * y)


def compute_error_convergence(N_list, analyzer, error_func):
        """
        Computa a lista de erros para diferentes ordens N usando a função de erro fornecida.

        Parâmetros:
        -----------
        N_list : iterable
            Lista de ordens polinomiais a testar.
        analyzer : objeto
            Objeto que contém os métodos de erro.
        error_func : callable
            Função de erro a ser chamada, como analyzer.reference_element_error ou analyzer.global_domain_error.

        Retorna:
        --------
        errors : list of float
            Lista com os erros computados.
        """
        errors = []
        for N in N_list:
            analyzer.N = N
            analyzer.N_plot = N + 1
            errors.append(error_func())

        df = analyzer.convergence_table(list(N_list), errors)
        print("\nConvergência no triângulo de referência:\n", df.to_string(index=False))
        return errors


def main():
    clear_terminal()
    mesh_data = mesh_rectangular_domain(PROBLEM, BOUNDARY, MATERIAL, h=1.0, view_mesh=False, mesh_info=False)
    mesh=Mesh2D(vx=mesh_data['VX'], vy=mesh_data['VY'], EToV=mesh_data['EToV'])

    analyzer = SpectralAnalyzer(mesh, n_order=3, exact_func=exact_u)
    plotter = SpectralVisualizer(analyzer, mesh)

    interpolate_data = analyzer.interpolate_solution()
    plotter.collocation_points()
    plotter.collocation_and_interpolated_points(interpolate_data)
    plotter.exact_solution_global(exact_u)
    plotter.interpolated_global_solution(analyzer)

    N_list = range(1, 15)
    # Para erro no elemento de referência
    ref_error = compute_error_convergence(N_list, analyzer, analyzer.reference_element_error)
    plotter.convergence_rate(N_list, ref_error, title="Convergência espectral no triângulo de referência")

    # Para erro global no domínio físico
    global_error = compute_error_convergence(N_list, analyzer, analyzer.global_domain_error)
    plotter.convergence_rate(N_list, global_error, title="Convergência espectral no domínio físico")


if __name__ == '__main__':
    main()
