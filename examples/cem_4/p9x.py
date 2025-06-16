import os
import sys

import scipy as sp

# Caminho do projeto
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..\..')))

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.tri import Triangulation
import matplotlib.ticker as ticker

from maxwell.dg.dg2d import *
from maxwell.driver import *
from maxwell.integrators.LSERK4 import *
from mesher.create_mesh import *
from maxwell.utils import *

BOUNDARY = [{'tag': 101, 'type': 'Dirichlet', 'value': None, 'name': 'circular_scatterer'}]

MATERIAL = [{'tag': 201, 'name': 'free_space', 'relative_magnetic_permeability': 1, 'relative_electric_permittivity': 1},
            {'tag': 301, 'name': 'PML_a','relative_magnetic_permeability': 1, 'relative_electric_permittivity': 1},
            {'tag': 302, 'name': 'PML_b','relative_magnetic_permeability': 1, 'relative_electric_permittivity': 1},
            {'tag': 303, 'name': 'PML_c','relative_magnetic_permeability': 1, 'relative_electric_permittivity': 1},
            {'tag': 304, 'name': 'PML_d','relative_magnetic_permeability': 1, 'relative_electric_permittivity': 1},
            {'tag': 401, 'name': 'PML_I','relative_magnetic_permeability': 1, 'relative_electric_permittivity': 1},
            {'tag': 402, 'name': 'PML_II','relative_magnetic_permeability': 1, 'relative_electric_permittivity': 1},
            {'tag': 403, 'name': 'PML_III','relative_magnetic_permeability': 1, 'relative_electric_permittivity': 1},
            {'tag': 404, 'name': 'PML_IV','relative_magnetic_permeability': 1, 'relative_electric_permittivity': 1}]

INFO_GRAPH = {'cell': True, 'nodes': True, 'edges': False, 'edges_numb': False, 'filepath': 'examplesData/inputs/cem_4p9/cem_4p9.svg'}

PROBLEM = {'name': 'cem_4p9', 'folder_name': 'cem_4p9',
            'flux_type': 'Centered',    # 'Upwind' or 'Centered'
            'cfl': 0.1,         # Número de Courant-Friedrichs-Lewy
            'bc': "PEC",        # Condição de contorno do problema: Perfect Electric Conductor (PEC)
            'n_order': 3,       # Ordem polinomial dos elementos
            'Lx': 2.0,          # Largura do domínio na direção x
            'Ly': 2.0,          # Largura do domínio na direção y
            'L': 1.0,           # Largura da camada de PML
            'x0': 1.0,          # Distância da fronteira interna da PML
            'pml_order': 3,     # Ordem polinomial do PML
            'R': 1E-4,          # Coeficiente de reflexão na interface do PML
}


class SpectralAnalyzer:
    def __init__(self, sp, dg_order):
        self.mesh = sp.mesh
        self.sp = sp
        self.N = dg_order                           # Ordem polinomial dos elementos
        self.N_plot = self.N                        # Ordem dos pontos de visualização

        # self.sigma_x, self.sigma_y = sp.build_sigma_fields(
        #     x_limit=PROBLEM['x0'],
        #     y_limit=PROBLEM['x0'],
        #     sigma0_x=-np.log(PROBLEM['R']),
        #     sigma0_y=-np.log(PROBLEM['R']),
        #     p=2
        # )


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


    def collocation_points(self, show_ids=False, sigma_data=False):
        """
        Plota os pontos de colocação no domínio físico com valores (σₓ, σᵧ) ao lado.
        Se show_ids for True, mostra os índices dos nós locais.
        """
        VX, VY, EToV = self.mesh.vx, self.mesh.vy, self.mesh.EToV
        K = EToV.shape[0]

        x_nodes, y_nodes = nodes_coordinates(self.N, self.mesh)
        markers = ['o', 's', '^', 'D', 'v', '>', '<', 'p', '*', 'h', 'X']
        markers = ['o']
        colors = plt.cm.tab10(np.linspace(0, 1, 10))

        for k in range(K):
            nodes = EToV[k]
            poly_x = VX[nodes]
            poly_y = VY[nodes]
            plt.plot(np.append(poly_x, poly_x[0]), np.append(poly_y, poly_y[0]), 'k-', lw=0.8)

        for k in range(K):
            xk = x_nodes[:, k]
            yk = y_nodes[:, k]
            mk = markers[k % len(markers)]
            ck = colors[k % len(colors)]
            plt.scatter(xk, yk, marker=mk, facecolors='none', color=ck, s=10, label=f'e{k}')

            # Centro do elemento
            va, vb, vc = EToV[k]
            xc = (VX[va] + VX[vb] + VX[vc]) / 3
            yc = (VY[va] + VY[vb] + VY[vc]) / 3

            for i in range(x_nodes.shape[0]):
                dx = xc - xk[i]
                dy = yc - yk[i]
                scale = 0.12
                xt = xk[i] + dx * scale
                yt = yk[i] + dy * scale

                if show_ids:
                    plt.text(xt, yt, str(i), fontsize=7, color='black', ha='center', va='center')

                if sigma_data:
                    label = f"({self.sigma_x[i, k]:.1f}, {self.sigma_y[i, k]:.1f})"
                    plt.text(xt, yt-0.06, label, fontsize=7, color='darkgreen', ha='center')

        plt.gca().set_aspect('equal')
        plt.xlabel("x")
        plt.ylabel("y")
        plt.title(fr"Pontos de colocação e campos $\sigma^x$, $\sigma^y$ ($N={self.N}$)")
        plt.grid(True, linestyle='--', linewidth=0.3)
        # plt.legend(fontsize=10, loc='center left', bbox_to_anchor=(1, 0.5))
        plt.tight_layout()
        plt.show()


    def plot_sigma_surfaces(self, field_labels=('σₓ(x, y)', 'σᵧ(x, y)')):
        """
        Plota as superfícies dos campos sigma_x e sigma_y no domínio computacional.

        Parâmetros:
        -----------
        sigma_x, sigma_y : ndarray
            Campos de amortecimento definidos nos pontos de colocação (shape: [Np, K]).
        field_labels : tuple
            Rótulos para os campos (eixos z).
        """
        from mpl_toolkits.mplot3d import Axes3D

        x_all, y_all, triangles_all = self.build_global_visualization_mesh()

        sigma_x_flat = self.sigma_x.ravel('F')
        sigma_y_flat = self.sigma_y.ravel('F')
        tri_global = Triangulation(x_all, y_all, triangles_all)

        fig = plt.figure(figsize=(14, 6))

        # σₓ plot
        ax1 = fig.add_subplot(1, 2, 1, projection='3d')
        surf1 = ax1.plot_trisurf(tri_global, sigma_x_flat, cmap='plasma', linewidth=0.2)
        ax1.set_title(f'Perfil de {field_labels[0]}')
        ax1.set_xlabel('x')
        ax1.set_ylabel('y')
        ax1.set_zlabel(field_labels[0])
        fig.colorbar(surf1, ax=ax1, shrink=0.6)

        # σᵧ plot
        ax2 = fig.add_subplot(1, 2, 2, projection='3d')
        surf2 = ax2.plot_trisurf(tri_global, sigma_y_flat, cmap='plasma', linewidth=0.2)
        ax2.set_title(f'Perfil de {field_labels[1]}')
        ax2.set_xlabel('x')
        ax2.set_ylabel('y')
        ax2.set_zlabel(field_labels[1])
        fig.colorbar(surf2, ax=ax2, shrink=0.6)

        for ax in (ax1, ax2):
            ax.view_init(elev=35, azim=-120)

        plt.tight_layout()
        plt.show()


    def set_gaussian_pulse(self, x, y, t, sigma=0.1, amplitude=1.0, c=1.0):
        """
        Pulso gaussiano esfericamente simétrico propagante no tempo.

        E_z(x, y, t) = A * exp(- (sqrt(x² + y²) - c*t)² / (2σ²))

        Parâmetros:
        - x, y: coordenadas espaciais
        - t: instante de tempo
        - sigma: largura do pulso
        - amplitude: valor máximo do campo
        - c: velocidade de propagação

        Retorna:
        - E_z(x, y, t): campo escalar
        """
        r = np.sqrt(x**2 + y**2)
        return amplitude * np.exp(-((r - c * t)**2) / (2 * sigma**2))


    def plot_analytical_time_evolution(self, t_list, N_grid=200):
        """
        Plota a evolução temporal da solução analítica do pulso gaussiano propagante
        sobre uma malha regular em [-1, 1]².

        Parâmetros:
        - t_list: lista de tempos em que o pulso será avaliado
        - N_grid: resolução da malha regular (default: 200)
        """
        # Malha regular em [-1, 1] x [-1, 1]
        x_lin = np.linspace(-1.0, 1.0, N_grid)
        y_lin = np.linspace(-1.0, 1.0, N_grid)
        X, Y = np.meshgrid(x_lin, y_lin)
        points = np.column_stack((X.ravel(), Y.ravel()))

        n_times = len(t_list)
        ncols = int(np.ceil(np.sqrt(n_times)))
        nrows = int(np.ceil(n_times / ncols))

        fig, axes = plt.subplots(nrows, ncols, figsize=(4.5 * ncols, 4 * nrows), constrained_layout=True)
        axes = axes.ravel()

        for i, t in enumerate(t_list):
            # Avaliação da solução propagante
            Ez_vals = self.set_gaussian_pulse(points[:, 0], points[:, 1], t)
            Ez_grid = Ez_vals.reshape(N_grid, N_grid)

            im = axes[i].imshow(Ez_grid, extent=[-1, 1, -1, 1], origin='lower',
                                cmap='viridis', aspect='equal', vmin=0, vmax=1)
            axes[i].set_title(f"t = {t:.1f} s")
            axes[i].set_xlabel("x")
            axes[i].set_ylabel("y")
            axes[i].set_xticks(np.linspace(-1, 1, 3))
            axes[i].set_yticks(np.linspace(-1, 1, 3))
            axes[i].grid(False)
            cbar = fig.colorbar(im, ax=axes[i])
            cbar.set_label(r'$E_z$ (V/m)', rotation=0, labelpad=12)

        # Desativar eixos extras
        for j in range(i + 1, len(axes)):
            axes[j].axis('off')

        fig.suptitle(r'$E_z(x,y,t)$ sobre domínio físico', fontsize=16)
        plt.show()


    def plot_global_solution(self, sp, field, field_label, ax=None):
        """
        Plota a solução DG diretamente nos pontos modais do método, com escala fixa.

        Parâmetros
        ----------
        sp : Maxwell2D
            Objeto com malha e campos do método DG.
        field : ndarray (Np x K)
            Solução numérica DG (Ez, Hx ou Hy).
        field_label : str
            Nome do campo a ser exibido no gráfico.
        ax : matplotlib.axes.Axes
            Subplot onde a solução será desenhada.
        """
        if ax is None:
            ax = plt.gca()

        x = sp.x.ravel(order='F')
        y = sp.y.ravel(order='F')
        uh = field.ravel(order='F')

        r, s = xy_to_rs(*set_nodes_in_equilateral_triangle(sp.n_order))
        tri_local = Triangulation(r, s).triangles

        Np, K = field.shape
        triangles_all = []
        for k in range(K):
            offset = k * Np
            tri_k = tri_local + offset
            triangles_all.extend(tri_k)
        tri_global = Triangulation(x, y, np.array(triangles_all))

        # Contorno 2D com escala de cores fixa
        tcf = ax.tricontourf(tri_global, uh, levels=100, cmap='viridis', vmin=0, vmax=1)
        cbar = plt.colorbar(tcf, ax=ax, label=r'$E_z$ (V/m)')
        cbar.formatter = ticker.FormatStrFormatter('%.1f')
        cbar.update_ticks()


        # Configuração visual
        ax.set_title(field_label)
        ax.set_xlabel("x")
        ax.set_ylabel("y")
        ax.set_xticks(np.linspace(-2, 2, 5))
        ax.set_yticks(np.linspace(-2, 2, 5))
        ax.set_aspect("equal")


    def run_with_plots(self, driver, T_final=1.5, n_frames=9):
        """
        Executa a simulação até T_final, coletando snapshots intermediários
        da variável Ez e plotando sua evolução temporal.

        Parâmetros:
        - driver: instância de MaxwellDriver
        - T_final: tempo final da simulação (s)
        - n_frames: número de quadros a exibir (default: 9)
        """
        dt = driver.dt
        sp = driver.sp
        t_plot_list = np.linspace(0, T_final, n_frames)

        Ez_history = []
        current_time = 0.0
        next_plot_index = 0
        t_next = t_plot_list[next_plot_index]

        while current_time < T_final and next_plot_index < len(t_plot_list):
            steps_to_next_plot = int(np.round((t_next - current_time) / dt))
            for _ in range(steps_to_next_plot):
                driver.step()
            current_time += steps_to_next_plot * dt

            # Armazena o campo Ez atual
            Ez_history.append(driver['Ez'].copy())
            next_plot_index += 1
            if next_plot_index < len(t_plot_list):
                t_next = t_plot_list[next_plot_index]

        # Subplots
        fig, axes = plt.subplots(3, 3, figsize=(16, 9), constrained_layout=True)
        axes = axes.ravel()

        for i, Ez in enumerate(Ez_history):
            self.plot_global_solution(sp, Ez, field_label=fr'$t={t_plot_list[i]:.1f}\,\mathrm{{s}}$', ax=axes[i])

        # Desativa eixos restantes
        for j in range(len(Ez_history), len(axes)):
            axes[j].axis('off')

        fig.suptitle(r" $E_z(x,y,t)$ with DG-FEM", fontsize=14)
        plt.show()



def main():
    clear_terminal()

    # Criar a malha retangular com PML
    mesh_data = mesh_rectangular_pml_domain(PROBLEM, BOUNDARY, MATERIAL, h=0.5, view_mesh=False)
    pml_mesh=Mesh2D(vx=mesh_data['VX'], vy=mesh_data['VY'], EToV=mesh_data['EToV'])
    # plot_triangular_mesh(INFO_GRAPH, mesh_data)

    # Definir a discretização espacial usando DG2D
    sp = Maxwell2D(n_order=PROBLEM['n_order'], mesh=pml_mesh, fluxType=PROBLEM['flux_type'])
    analyzer = SpectralAnalyzer(sp, dg_order=PROBLEM['n_order'])

    # Plotar pontos de colocação com os valores de σₓ e σᵧ
    # analyzer.collocation_points(show_ids=False, sigma_data=False)

    # Visualizar os campos σₓ e σᵧ como superfícies 3D
    # analyzer.plot_sigma_surfaces()

    # Solução analítica no tempo
    # t_list = np.linspace(0, 1.5, 9)
    # analyzer.plot_analytical_time_evolution(t_list)
    
    # Initialize the solver
    driver = MaxwellDriver(sp, CFL=PROBLEM['cfl'])
    driver['Ez'][:] = analyzer.set_gaussian_pulse(sp.x, sp.y, 0)

    # Simulação 
    analyzer.run_with_plots(driver, T_final=1.5, n_frames=9)
    

if __name__ == '__main__':
    main()
