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

BOUNDARY = [{'tag': 101, 'type': 'Dirichlet', 'value': 0.0, 'name': 'Ez_0'}]
MATERIAL = [{'tag': 201, 'name': 'free_space', 'relative_magnetic_permeability': 1, 'relative_electric_permittivity': 1}]

def exact_u(x, y):
    return np.cos(np.pi * x) * np.cos(np.pi * y)


class SpectralAnalyzer:
    def __init__(self, problem_data, mesh_data, exact_func):
        self.mesh_data = mesh_data
        self.exact_func = exact_func
        self.N = problem_data['n_order']

    def build_reference_interpolation(self):
        """
        Gera os nós de colocação no triângulo de referência e constrói a matriz de Vandermonde  
        para interpolação de uma função exata definida no triângulo de referência.
        """
        N = self.N
        x, y = set_nodes_in_equilateral_triangle(N)
        r, s = xy_to_rs(x, y)
        V = vandermonde(N, r, s)
        u_nodes = self.exact_func(r, s)
        coeff = np.linalg.inv(V) @ u_nodes
        return coeff, V
    
    def interpolate_solution_on_physical_domain(self):
        """
        Interpola a solução exata u(x,y) para pontos físicos de visualização
        distribuídos em todos os elementos da malha.

        Retorna:
            dict contendo 'x_plot', 'y_plot', 'u_plot' com shape (Np_plot, K)
        """
        N = self.N
        mesh = self.mesh_data
        VX, VY, EToV = mesh['VX'], mesh['VY'], mesh['EToV']
        K = EToV.shape[0]

        # Pontos de visualização no triângulo de referência
        x_ref, y_ref = set_nodes_in_equilateral_triangle(N + 2)
        r_plot, s_plot = xy_to_rs(x_ref, y_ref)
        V_plot = vandermonde(N, r_plot, s_plot)
        Np_plot = len(r_plot)

        x_plot = np.zeros((Np_plot, K))
        y_plot = np.zeros((Np_plot, K))
        u_plot = np.zeros((Np_plot, K))

        for k in range(K):
            va, vb, vc = EToV[k]
            xk = 0.5 * (-(r_plot + s_plot) * VX[va] + (1 + r_plot) * VX[vb] + (1 + s_plot) * VX[vc])
            yk = 0.5 * (-(r_plot + s_plot) * VY[va] + (1 + r_plot) * VY[vb] + (1 + s_plot) * VY[vc])
            x_plot[:, k] = xk
            y_plot[:, k] = yk

            # Avalia u nos pontos nodais do elemento
            r_nodes, s_nodes = xy_to_rs(*set_nodes_in_equilateral_triangle(N))
            x_nodes = 0.5 * (-(r_nodes + s_nodes) * VX[va] + (1 + r_nodes) * VX[vb] + (1 + s_nodes) * VX[vc])
            y_nodes = 0.5 * (-(r_nodes + s_nodes) * VY[va] + (1 + r_nodes) * VY[vb] + (1 + s_nodes) * VY[vc])
            u_nodes = self.exact_func(x_nodes, y_nodes)

            V = vandermonde(N, r_nodes, s_nodes)
            coeff_k = np.linalg.solve(V, u_nodes)

            u_plot[:, k] = V_plot @ coeff_k

        return {
            'x_plot': x_plot,
            'y_plot': y_plot,
            'u_plot': u_plot,
            'tri_ref': Triangulation(r_plot, s_plot)
        }

    def compute_L2_error(self, u_approx, wq):
        return np.sqrt(np.sum(wq * (self.exact_func - u_approx)**2))

    def compute_error_reference(self):
        N = self.N
        coeff, V = self.build_reference_interpolation(N)
        Nq = N + 3
        xq, yq = set_nodes_in_equilateral_triangle(Nq)
        rq, sq = xy_to_rs(xq, yq)
        Vq = vandermonde(N, rq, sq)
        u_interp = Vq @ coeff
        u_exact = self.exact_func(rq, sq)
        wq = np.full_like(rq, 2.0 / len(rq))
        return self.compute_L2_error(u_exact, u_interp, wq)

    def compute_error_global(self):
        N = self.N
        VX, VY, EToV = self.mesh_data['VX'], self.mesh_data['VY'], self.mesh_data['EToV']
        coeff, V = self.build_reference_interpolation(N)
        Nq = N + 3
        xq, yq = set_nodes_in_equilateral_triangle(Nq)
        rq, sq = xy_to_rs(xq, yq)
        Vq = vandermonde(N, rq, sq)
        u_interp_local = Vq @ coeff
        wq = np.full_like(rq, 2.0 / len(rq))
        error_total = 0.0
        for k in range(EToV.shape[0]):
            va, vb, vc = EToV[k]
            xk = 0.5 * (-(rq + sq) * VX[va] + (1 + rq) * VX[vb] + (1 + sq) * VX[vc])
            yk = 0.5 * (-(rq + sq) * VY[va] + (1 + rq) * VY[vb] + (1 + sq) * VY[vc])
            Jk = abs((VX[vb] - VX[va]) * (VY[vc] - VY[va]) - (VX[vc] - VX[va]) * (VY[vb] - VY[va])) / 2
            u_exact_k = self.exact_func(xk, yk)
            error_total += np.sum(wq * (u_exact_k - u_interp_local)**2) * Jk
        return np.sqrt(error_total)

    def convergence_table(self, N_list, errors):
        rates = [np.nan]
        for i in range(1, len(N_list)):
            e0, e1 = errors[i - 1], errors[i]
            rate = np.log(e0 / e1) / np.log(N_list[i] / N_list[i - 1]) if e0 > 0 and e1 > 0 else np.nan
            rates.append(round(rate, 2))
        return pd.DataFrame({"N": N_list, "Erro L2": errors, "Taxa": rates})


class SpectralVisualizer():
    def __init__(self, analyzer, mesh_data):
        self.mesh_data = mesh_data
        self.N = analyzer.N

    def build_global_visualization_mesh(self, N_plot):
        """
        Gera uma malha global de visualização a partir da malha física, com conectividades contínuas.

        Retorna:
        --------
        x_all, y_all : ndarray
            Todos os pontos físicos (flattened) da visualização.
        triangles_all : ndarray
            Conectividade global dos triângulos (para Triangulation).
        """
        VX, VY, EToV = self.mesh_data['VX'], self.mesh_data['VY'], self.mesh_data['EToV']
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

    def plot_collocation_points(self, show_ids=False):
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
        VX = self.mesh_data['VX']
        VY = self.mesh_data['VY']
        EToV = self.mesh_data['EToV']
        K = EToV.shape[0]

        x_nodes, y_nodes = nodes_coordinates_from_gmsh(self.N, self.mesh_data)
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

    def plot_collocation_and_interpolated_points(self, interpolate_data):
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
        x_nodes, y_nodes = nodes_coordinates_from_gmsh(self.N, self.mesh_data)
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
        VX = self.mesh_data['VX']
        VY = self.mesh_data['VY']
        EToV = self.mesh_data['EToV']
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

    def plot_exact_solution_global(self, exact_func, field_label='u'):
        """
        Plota a solução exata u(x,y) sobre uma malha global contínua.
        """
        N_plot = max(self.N + 3, 6)  # Número de pontos para visualização
        x_all, y_all, triangles_all = self.build_global_visualization_mesh(N_plot)
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

    def plot_interpolated_solution_global_continuous(self, analyzer, field_label='u'):
        """
        Plota a solução interpolada u(x,y) globalmente com continuidade entre elementos.
        Utiliza avaliação local por elemento e unifica a malha e a solução interpolada.
        """
        N = self.N
        VX, VY, EToV = self.mesh_data['VX'], self.mesh_data['VY'], self.mesh_data['EToV']
        K = EToV.shape[0]

        # Pontos de visualização no triângulo de referência
        x_ref, y_ref = set_nodes_in_equilateral_triangle(N + 2)
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
    


def main():
    clear_terminal()
    PROBLEM = {'name': 'cem_4p7_u(x)', 'folder_name': 'cem_4p7_u(x)', 'L': 2.0, 'n_order': 4}
    mesh_data = mesh_rectangular_domain(PROBLEM, BOUNDARY, MATERIAL, h=1.0, view_mesh=False, mesh_info=False)

    analyzer = SpectralAnalyzer(PROBLEM, mesh_data, exact_func=exact_u)
    visualizer = SpectralVisualizer(analyzer, mesh_data)

    interpolate_data = analyzer.interpolate_solution_on_physical_domain()
    visualizer.plot_collocation_points()
    visualizer.plot_collocation_and_interpolated_points(interpolate_data)
    visualizer.plot_exact_solution_global(exact_u)
    visualizer.plot_interpolated_solution_global_continuous(analyzer)

    # N_list = range(1, 15)
    # errors = [analyzer.compute_error_reference(N) for N in N_list]
    # df = analyzer.convergence_table(list(N_list), errors)
    # print("\nConvergência no triângulo de referência:\n", df.to_string(index=False))
    # visualizer.plot_convergence_curve(N_list, errors, title="Convergência espectral no triângulo de referência")


if __name__ == '__main__':
    main()
