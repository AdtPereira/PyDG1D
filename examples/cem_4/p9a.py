#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# pylint: disable=E0611, C0413, C0103, W0401, W0614, W0622

import os
import sys

# Adiciona a raiz do projeto ao PYTHONPATH
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..\..')))

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.tri import Triangulation
from scipy.interpolate import griddata

from maxwell.dg.dg2d_pml import *
from maxwell.driver import *
from maxwell.integrators.LSERK4 import *
from mesher.create_mesh import *
from maxwell.utils import *

BOUNDARY = [{'tag': 101, 'type': 'Dirichlet', 'value': None, 'name': 'outermost_domain'}]
PROBLEM = {
    'DIM': 2,
    'name': 'p9',
    'folder_name': 'cem_4',
    'description': 'Análise espectral de um pulso gaussiano propagante em 2D',
    'flux_type': 'Upwind',    # 'Upwind' or 'Centered'
    'bc': "PEC",        # Condição de contorno do problema: Perfect Electric Conductor (PEC)
    'cfl': 0.1,         # Número de Courant-Friedrichs-Lewy
    'n_order': 3,       # Ordem polinomial dos elementos
    't_final': 5.0      # Tempo final da simulação
}


class SpectralAnalyzer:
    def __init__(self, problem, mesh):
        self.N = problem['n_order']
        self.L = 2*problem['Lx']
        self.t_final = problem['t_final']
        self.mesh = mesh


    def gaussian_pulse(self, x, y, t, A0=1.0, sigma_r=0.125, sigma_t=1.0, t0=0.0):
        """
        Calcula a amplitude de um pulso gaussiano 2D no domínio do tempo.

        O pulso é centrado espacialmente em (0,0) e tem seu pico temporal em t0.

        Args:
            x (float or np.ndarray): Coordenada(s) x.
            y (float or np.ndarray): Coordenada(s) y.
            t (float): Instante de tempo.
            A0 (float, optional): Amplitude máxima do pulso. Default é 1.0.
            sigma_r (float, optional): Desvio padrão espacial (largura do pulso).
                                    Default é 0.1.
            sigma_t (float, optional): Desvio padrão temporal (duração do pulso).
                                    Default é 1.0.
            t0 (float, optional): Tempo de pico do pulso. Default é 0.0.

        Returns:
            float or np.ndarray: A amplitude do pulso no(s) ponto(s) (x,y) no instante t.
        """
        termo_espacial = np.exp(-(x**2 + y**2) / (2 * sigma_r**2))
        termo_temporal = np.exp(-((t - t0)**2) / (2 * sigma_t**2))
        return A0 * termo_espacial * termo_temporal


    def gaussian_pulse_derivative(self, x, y, t, sigma=0.2, amplitude=1.0, c=1.0):
        """
        Derivada temporal do pulso gaussiano.

        dE_z/dt = A * [c*(r - c*t)/σ²] * exp(- (r - c*t)² / (2σ²))

        Parâmetros:
        - x, y: coordenadas espaciais
        - t: instante de tempo
        - sigma: largura do pulso
        - amplitude: valor máximo do campo original
        - c: velocidade de propagação

        Retorna:
        - dE_z/dt: valor da derivada temporal do campo
        """
        r = np.sqrt(x**2 + y**2)
        
        # Termo que surge da regra da cadeia na diferenciação
        pre_factor = amplitude * (c * (r - c * t)) / (sigma**2)
        
        # Parte exponencial do pulso gaussiano original
        gaussian_part = np.exp(-((r - c * t)**2) / (2 * sigma**2))
        
        return pre_factor * gaussian_part


    def single_test(self, driver):
        """
        Teste simples da classe SpectralAnalyzer.
        """
        # Set the time integrator
        driver.run_until(self.t_final)
        n_steps = int(np.ceil(self.t_final / driver.dt))
        print(f"\n🌐 Passo de tempo da simulação: {(1E3*driver.dt):.2f} ms")
        print(f"🌐 Número de passos de tempo: {n_steps}")
        print(f"🌐 Tempo final da simulação: {driver.timeIntegrator.time:.3f} s")

        ez_expected = self.gaussian_pulse(driver.sp.x, driver.sp.y, driver.timeIntegrator.time)
        R = np.corrcoef(ez_expected, driver['Ez'])
        print(f"🌐 Correlação entre solução esperada e solução numérica: {R[0,1]:.2f}")

        # Plota a solução numérica Ez diretamente nos pontos modais e sobre malha regular
        # self.plot_modal_global_solution(driver.sp, driver['Ez'], field_label='E_z', show_mesh=False)
        self.plot_griddata_global_solution(driver.sp, driver['Ez'], field_label='E_z', N_grid=300)


    def plot_collocation_points(self, show_ids=False, sigma_data=False):
        """
        Plota os pontos de colocação no domínio físico com valores (σₓ, σᵧ) ao lado.
        Se show_ids for True, mostra os índices dos nós locais.
        """
        VX, VY, EToV = self.mesh.vx, self.mesh.vy, self.mesh.EToV
        K = EToV.shape[0]

        x_nodes, y_nodes = nodes_coordinates(self.N, self.mesh)
        markers = ['o', 's', '^', 'D', 'v', '>', '<', 'p', '*', 'h', 'X']
        markers = ['o']
        plt.figure(figsize=(10, 8))
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
            plt.scatter(xk, yk, marker=mk, facecolors='none', color=ck, s=4, label=f'e{k}')

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
        plt.title(fr"Pontos de colocação for $N={self.N}$")
        plt.grid(True, linestyle='--', linewidth=0.3)
        if K < 10:
            plt.legend(fontsize=10, loc='center left', bbox_to_anchor=(1, 0.5))
        plt.tight_layout()


    def plot_analytical_field_at_final_time(self, Nx=200, Ny=200):
        """
        Plota a solução analítica Ez sobre o domínio físico [-π, π]^2
        em t = 3T, com dois subplots: mapa de cores 2D e superfície 3D.
        """
        # Domínio
        L = self.L
        x = np.linspace(-L/2, L/2, Nx)
        y = np.linspace(-L/2, L/2, Ny)
        X, Y = np.meshgrid(x, y)

        # Avaliação da solução
        Ez = self.gaussian_pulse(X, Y, self.t_final)

        # Criação dos subplots
        fig = plt.figure(figsize=(12, 5))
        fig.suptitle(rf'Solução Analítica $E_z(x, y, t = {self.t_final:.1f} s)$', fontsize=14)

        # Subplot 1: Heatmap 2D
        ax1 = fig.add_subplot(1, 2, 1)
        pcm = ax1.pcolormesh(X, Y, Ez, shading='auto', cmap='viridis')
        ax1.set_xlabel('x')
        ax1.set_ylabel('y')
        ax1.axis('equal')
        fig.colorbar(pcm, ax=ax1, label=r'$E_z$')

        # Subplot 2: Superfície 3D
        ax2 = fig.add_subplot(1, 2, 2, projection='3d')
        surf = ax2.plot_surface(X, Y, Ez, cmap='viridis', edgecolor='none')
        ax2.set_xlabel('x')
        ax2.set_ylabel('y')
        ax2.set_zlabel(r'$E_z$')
        fig.colorbar(surf, ax=ax2, shrink=0.6, aspect=12)

        plt.tight_layout()
        plt.show()


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
            Ez_vals = self.gaussian_pulse(points[:, 0], points[:, 1], t)
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


    def plot_modal_global_solution(self, sp, field, field_label='u', show_mesh=True):
        """
        Plota a solução DG diretamente nos pontos modais do método, sem interpolação.

        Parâmetros
        ----------
        sp : Maxwell2D
            Objeto com malha e campos do método DG.
        field : ndarray (Np x K)
            Solução numérica DG (Ez, Hx ou Hy).
        field_label : str
            Nome do campo a ser exibido no gráfico.
        show_mesh : bool
            Se True, sobrepõe as bordas dos elementos triangulares.
        """
        # Reorganiza os dados em arrays planos
        x = sp.x.ravel(order='F')
        y = sp.y.ravel(order='F')
        uh = field.ravel(order='F')

        r, s = xy_to_rs(*set_nodes_in_equilateral_triangle(sp.n_order))
        tri_local = Triangulation(r, s).triangles

        # Constrói triangulação global
        Np, K = field.shape
        triangles_all = []
        for k in range(K):
            offset = k * Np
            tri_k = tri_local + offset
            triangles_all.extend(tri_k)
        tri_global = Triangulation(x, y, np.array(triangles_all))

        # --- Gráficos ---
        fig = plt.figure(figsize=(14, 6))
        ax1 = fig.add_subplot(1, 2, 1)
        ax2 = fig.add_subplot(1, 2, 2, projection='3d')

        # Contorno 2D
        tcf = ax1.tricontourf(tri_global, uh, levels=100, cmap='viridis')
        fig.colorbar(tcf, ax=ax1, label=fr'{field_label}(x,y)')
        ax1.set_title(fr'Solução DG ${field_label}(x, y)$ nos pontos modais')
        ax1.set_xlabel("x")
        ax1.set_ylabel("y")
        ax1.set_aspect("equal")

        if show_mesh:
            VX, VY, EToV = sp.mesh.vx, sp.mesh.vy, sp.mesh.EToV
            for tri in EToV:
                px = [VX[tri[0]], VX[tri[1]], VX[tri[2]], VX[tri[0]]]
                py = [VY[tri[0]], VY[tri[1]], VY[tri[2]], VY[tri[0]]]
                ax1.plot(px, py, color='k', linewidth=0.5)

        # Superfície 3D
        ax2.plot_trisurf(tri_global, uh, cmap='viridis', linewidth=0.2)
        ax2.set_title(fr'Solução DG ${field_label}(x, y)$ - Superfície')
        ax2.set_xlabel("x")
        ax2.set_ylabel("y")
        ax2.set_zlabel(fr'{field_label}(x, y)')
        ax2.view_init(elev=35, azim=-120)
        plt.tight_layout()
        plt.show()


    def plot_griddata_global_solution(self, sp, field, field_label='u', N_grid=300):
        """
        Interpola a solução DG usando griddata sobre uma malha regular e plota
        mapa de cores e superfície 3D.

        Parâmetros:
        -----------
        sp : Maxwell2D
            Objeto do método DG contendo a malha e coordenadas dos pontos.
        field : ndarray (Np x K)
            Solução numérica nos pontos modais do DG.
        field_label : str
            Nome do campo para rótulo do gráfico.
        N_grid : int
            Número de pontos em cada direção da malha regular.
        """
        # Dados da solução nos pontos modais
        x = sp.x.ravel(order='F')
        y = sp.y.ravel(order='F')
        u = field.ravel(order='F')

        # Domínio de interpolação regular
        x_min, x_max = x.min(), x.max()
        y_min, y_max = y.min(), y.max()
        xi = np.linspace(x_min, x_max, N_grid)
        yi = np.linspace(y_min, y_max, N_grid)
        XI, YI = np.meshgrid(xi, yi)

        # Interpolação
        UI = griddata(points=(x, y), values=u, xi=(XI, YI), method='cubic')

        # Gráficos
        fig = plt.figure(figsize=(14, 6))
        fig.suptitle(fr'Interpolação com `griddata` da solução ${field_label}(x, y)$', fontsize=14)

        # Subplot 1: Mapa de cores
        ax1 = fig.add_subplot(1, 2, 1)
        img = ax1.pcolormesh(XI, YI, UI, shading='auto', cmap='viridis')
        fig.colorbar(img, ax=ax1, label=field_label)
        ax1.set_title("Mapa de Cores")
        ax1.set_xlabel("x")
        ax1.set_ylabel("y")
        ax1.set_aspect("equal")

        # Subplot 2: Superfície 3D
        ax2 = fig.add_subplot(1, 2, 2, projection='3d')
        surf = ax2.plot_surface(XI, YI, UI, cmap='viridis', edgecolor='none')
        ax2.set_title("Superfície 3D")
        ax2.set_xlabel("x")
        ax2.set_ylabel("y")
        ax2.set_zlabel(field_label)
        ax2.view_init(elev=35, azim=-120)

        plt.tight_layout()
        plt.show()
    

    def plot_numerical_time_evolution(self, times_snap, field_snapshots, driver, n_snaps=9, N_grid=300):
        """
        Plota a evolução temporal da solução numérica DG interpolada via griddata.

        Parâmetros:
        -----------
        field_snapshots : list of ndarray
            Lista com soluções DG (shape: Np x K) para cada instante de tempo.
        sp : Maxwell2D
            Objeto DG contendo malha e pontos modais.
        t_list : list of float
            Lista dos tempos correspondentes a cada solução.
        N_grid : int
            Número de pontos por eixo na malha regular para interpolação.
        """
        
        x, y = driver.sp.x.ravel(order='F'), driver.sp.y.ravel(order='F')
        xi = np.linspace(x.min(), x.max(), N_grid)
        yi = np.linspace(y.min(), y.max(), N_grid)
        XI, YI = np.meshgrid(xi, yi)

        ncols = int(np.ceil(np.sqrt(n_snaps)))
        nrows = int(np.ceil(n_snaps / ncols))

        fig, axes = plt.subplots(nrows, ncols, figsize=(4.5 * ncols, 4 * nrows), constrained_layout=True)
        axes = axes.ravel()

        vmin = min(np.min(f) for f in field_snapshots)
        vmax = max(np.max(f) for f in field_snapshots)

        for i, (Ez_k, t_k) in enumerate(zip(field_snapshots, times_snap)):
            u_k = Ez_k.ravel(order='F')
            Ui = griddata(points=(x, y), values=u_k, xi=(XI, YI), method='cubic')
            im = axes[i].imshow(
                Ui,
                extent=[x.min(), x.max(), y.min(), y.max()],
                origin='lower', cmap='gray_r', aspect='equal', vmin=vmin, vmax=vmax)

            axes[i].set_title(f"t = {t_k:.2f} s")
            axes[i].set_xlabel("x")
            axes[i].set_ylabel("y")
            axes[i].set_xticks(np.linspace(x.min(), x.max(), 5))
            axes[i].set_yticks(np.linspace(y.min(), y.max(), 5))
            axes[i].grid(False)
            cbar = fig.colorbar(im, ax=axes[i])
            cbar.set_label(r'$E_z$ (V/m)', rotation=90, labelpad=10)

        for j in range(i + 1, len(axes)):
            axes[j].axis('off')

        fig.suptitle(r'DG-FEM Solution for $E_z(x, y, t)$', fontsize=12)


    def Plot_ez_transient_analysis(self, times, ez_values):
        """
        Plota a evolução temporal do campo elétrico Ez em um ponto específico.

        Parâmetros:
        -----------
        times : list of float
            Lista de tempos correspondentes aos valores de Ez.
        ez_values : list of float
            Lista de valores do campo elétrico Ez no ponto monitorado.
        """
        plt.figure(figsize=(10, 5))
        plt.plot(times, ez_values, label=r'$E_z(x=1, y=0, t)$', color='blue')
        plt.xlabel("Tempo (s)")
        plt.ylabel(r'$E_z$ (V/m)')
        plt.title(r'Evolução do Campo Elétrico $E_z$ no Ponto (x=1, y=0)')
        plt.grid(True, linestyle='--', alpha=0.7)
        plt.legend()
        plt.tight_layout()


    def plot_sigma_surfaces(self, sp, N_grid=300):
        """
        Plota as superfícies dos campos sigma_x e sigma_y no domínio computacional.

        Parâmetros:
        -----------
        sigma_x, sigma_y : ndarray
            Campos de amortecimento definidos nos pontos de colocação (shape: [Np, K]).
        field_labels : tuple
            Rótulos para os campos (eixos z).
        """
        
        # Dados da solução nos pontos modais
        x = sp.x.ravel(order='F')
        y = sp.y.ravel(order='F')
        sigma_x = sp.sigma_x.ravel(order='F')
        sigma_y = sp.sigma_y.ravel(order='F')
        dsigma_dx = sp.dsigma_dx.ravel(order='F')
        dsigma_dy = sp.dsigma_dy.ravel(order='F')

        # Domínio de interpolação regular
        x_min, x_max = x.min(), x.max()
        y_min, y_max = y.min(), y.max()
        xi = np.linspace(x_min, x_max, N_grid)
        yi = np.linspace(y_min, y_max, N_grid)
        XI, YI = np.meshgrid(xi, yi)

        # Interpolação
        UIx = griddata(points=(x, y), values=sigma_x, xi=(XI, YI), method='cubic')
        UIy = griddata(points=(x, y), values=sigma_y, xi=(XI, YI), method='cubic')
        dUIx = griddata(points=(x, y), values=dsigma_dx, xi=(XI, YI), method='cubic')
        dUIy = griddata(points=(x, y), values=dsigma_dy, xi=(XI, YI), method='cubic')

        # Gráficos
        fig = plt.figure(figsize=(14, 6))
        fig.suptitle(fr'Interpolação do peril de condutividade da camada PML', fontsize=12)

        # Subplot 1: σₓ plot
        ax1 = fig.add_subplot(2, 2, 1, projection='3d')
        surf = ax1.plot_surface(XI, YI, UIx, cmap='viridis', edgecolor='none', vmin=0, vmax=10)
        ax1.set_title(r"$\sigma_x(x, y)$")
        ax1.set_xlabel("x")
        ax1.set_ylabel("y")
        ax1.set_zlabel(r"$\sigma_x$ (S/m)")
        fig.colorbar(surf, ax=ax1, shrink=0.6, aspect=12)
        #ax1.view_init(elev=35, azim=-120)

        # Subplot 2: σᵧ plot
        ax2 = fig.add_subplot(2, 2, 2, projection='3d')
        surf = ax2.plot_surface(XI, YI, UIy, cmap='viridis', edgecolor='none', vmin=0, vmax=10)
        ax2.set_title(r"$\sigma_y(x, y)$")
        ax2.set_xlabel("x")
        ax2.set_ylabel("y")
        ax2.set_zlabel(r"$\sigma_y$ (S/m)")
        fig.colorbar(surf, ax=ax2, shrink=0.6, aspect=12)
        #ax2.view_init(elev=35, azim=-120)

        # Subplot 3: dσₓ/dx plot
        ax3 = fig.add_subplot(2, 2, 3, projection='3d') 
        surf = ax3.plot_surface(XI, YI, dUIx, cmap='viridis', edgecolor='none', vmin=0, vmax=10)
        ax3.set_title(r"$\frac{\partial }{\partial x} \sigma_x (x, y)$")  
        ax3.set_xlabel("x")
        ax3.set_ylabel("y")
        ax3.set_zlabel(r"$\frac{\partial \sigma_x}{\partial x}$ (S/m²)")
        fig.colorbar(surf, ax=ax3, shrink=0.6, aspect=12)
        #ax3.view_init(elev=35, azim=-120)

        # Subplot 4: dσᵧ/dy plot
        ax4 = fig.add_subplot(2, 2, 4, projection='3d')
        surf = ax4.plot_surface(XI, YI, dUIy, cmap='viridis', edgecolor='none', vmin=0, vmax=10)
        ax4.set_title(r"$\frac{\partial }{\partial y} \sigma_y (x, y)$")
        ax4.set_xlabel("x")
        ax4.set_ylabel("y")
        ax4.set_zlabel(r"$\frac{\partial \sigma_y}{\partial y}$ (S/m²)")
        fig.colorbar(surf, ax=ax4, shrink=0.6, aspect=12)
        #ax4.view_init(elev=35, azim=-120)

        plt.tight_layout()
        plt.show()


def maxwell_at_free_space(problem, Lx):
    """
    Monitors and plots the electric field Ez at the PML interface (x=1, y=0) over time.
    """
    print(f"\n🌐 Inicializando o analisador espectral para {problem['name']} no espaço livre ...")
    # Define o domínio físico
    material = [{'tag': 201, 'name': 'free_space', 'relative_magnetic_permeability': 1, 'relative_electric_permittivity': 1}] 
    problem['Lx'], problem['Ly'] = Lx, Lx
    mesh_data = mesh_rectangular_domain(problem, BOUNDARY, material, h=0.5, view_mesh=False)
    mesh_domain = Mesh2D(vx=mesh_data['VX'], vy=mesh_data['VY'], EToV=mesh_data['EToV'], boundary_label=problem['bc'])
    
    # Initialize the solver
    sp = Maxwell2D(n_order=problem['n_order'], mesh=mesh_domain, fluxType=problem['flux_type'])

    return mesh_domain, sp


def maxwell_at_pml_space(problem, L):
    """
    Monitors and plots the electric field Ez at the PML interface (x=1, y=0) over time.
    """
    print(f"\n🌐 Inicializando o analisador espectral para {problem['name']} truncado com PML ...")
    # Define o domínio físico
    material = [{'tag': 201, 'name': 'free_space', 'relative_magnetic_permeability': 1, 'relative_electric_permittivity': 1},
                {'tag': 301, 'name': 'PML_a','relative_magnetic_permeability': 1, 'relative_electric_permittivity': 1},
                {'tag': 302, 'name': 'PML_b','relative_magnetic_permeability': 1, 'relative_electric_permittivity': 1},
                {'tag': 303, 'name': 'PML_c','relative_magnetic_permeability': 1, 'relative_electric_permittivity': 1},
                {'tag': 304, 'name': 'PML_d','relative_magnetic_permeability': 1, 'relative_electric_permittivity': 1},
                {'tag': 401, 'name': 'PML_I','relative_magnetic_permeability': 1, 'relative_electric_permittivity': 1},
                {'tag': 402, 'name': 'PML_II','relative_magnetic_permeability': 1, 'relative_electric_permittivity': 1},
                {'tag': 403, 'name': 'PML_III','relative_magnetic_permeability': 1, 'relative_electric_permittivity': 1},
                {'tag': 404, 'name': 'PML_IV','relative_magnetic_permeability': 1, 'relative_electric_permittivity': 1}]
    
    pml_project = {
        'x0': 1.0,          # Fronteira direita da PML
        'L': 1.0,           # Largura da camada da PML
        'pml_order': 2,     # Ordem polinomial da PML
        'R': 1E-4,          # Coeficiente de reflexão na interface da PML
    }

    pml_project['L'] = L
    
    mesh_data = mesh_rectangular_pml_domain(problem, pml_project, BOUNDARY, material, h=0.3, view_mesh=False)
    
    mesh_domain = Mesh2D(vx=mesh_data['VX'], vy=mesh_data['VY'], EToV=mesh_data['EToV'], boundary_label=problem['bc'])
    
    # Initialize the solver
    sp = Maxwell2D(n_order=problem['n_order'], mesh=mesh_domain, fluxType=problem['flux_type'], pml_design=pml_project)

    return mesh_domain, sp


def transient_analysis(sp, mesh_domain, point=(1, 0), n_samples=500, n_snaps=9):
    """
    Monitors and plots the electric field Ez at the PML interface (x=1, y=0) over time.
    """
    print(f"\n🌐 Inicializando o analisador espectral para {PROBLEM['name']} ...")
    sa = SpectralAnalyzer(PROBLEM, mesh_domain)
    # sa.plot_collocation_points(show_ids=False, sigma_data=False)

    # 1. Inicializar o campo elétrico Ez 
    driver = MaxwellDriver(sp, CFL=PROBLEM['cfl'])
    driver['Ez'][:] = sa.gaussian_pulse(sp.x, sp.y, 0)
    #driver['Ez'][:] = sa.gaussian_pulse_derivative(sp.x, sp.y, 0)

    # # 2. Calcular o raio 'r' para cada nó da malha
    # #    (adicionar um valor pequeno, epsilon, para evitar divisão por zero na origem)
    # epsilon = 1e-12
    # r = np.sqrt(sp.x**2 + sp.y**2) + epsilon

    # # 3. Calcular Hx e Hy usando a relação correta
    # #    Neste caso, a impedância eta = 1.0
    # eta = 1.0
    # Hx = (sp.y / r) * (Ez / eta)
    # Hy = (-sp.x / r) * (Ez / eta)

    # # 4. Atribuir os valores calculados aos campos do solver
    # driver['Hx'][:] = Hx
    # driver['Hy'][:] = Hy

    # Find the node closest to (x=1, y=0)
    x_point, y_point = point
    distances = np.sqrt((sp.x - x_point)**2 + (sp.y - y_point)**2)
    node_idx, elem_idx = np.unravel_index(np.argmin(distances), distances.shape)

    print(f"\n🌐 Monitorando o campo elétrico Ez no ponto ({point[0]}, {point[1]}) ...")
    print(f"  Element index: {elem_idx}, Node index: {node_idx}")
    print(f"  Actual coordinates: ({sp.x[node_idx, elem_idx]:.4f}, {sp.y[node_idx, elem_idx]:.4f})")

    times, ez_values, ez_snap, times_snap = [], [], [], []
    n_steps = int(sa.t_final / driver.dt)

    # Condição para capturar os dados da série temporal (ez_values)
    sample_every = max(1, n_steps // n_samples)

    # Determina os passos exatos em que os snapshots do campo serão capturados.
    if n_snaps > 0:
        # Cria um conjunto (set) com os índices dos passos para uma busca rápida (O(1)) no loop.
        snapshot_steps = set(np.round(np.linspace(0, n_steps - 1, n_snaps)).astype(int))
    else:
        snapshot_steps = set() # Conjunto vazio se nenhum snapshot for necessário.

    print(f"\n🌐 Número de passos de tempo: {n_steps}, amostrando a cada {sample_every} passos. {n_samples} amostras no total")
        
    for step in range(n_steps):
        driver.step()
        
        # Condição para capturar a série temporal em um ponto.
        if step % sample_every == 0:
            times.append(driver.timeIntegrator.time)
            ez_values.append(driver['Ez'][node_idx, elem_idx])

        # Condição para capturar o snapshot do campo completo.
        if step in snapshot_steps:
            times_snap.append(driver.timeIntegrator.time)
            ez_snap.append(driver['Ez'].copy())

    print(f"\n🌐 Tempo final da simulação: {driver.timeIntegrator.time:.3f} s")
    print(f"🌐 Passo de tempo da simulação: {(1E3*driver.dt):.2f} ms")
    print(f"🌐 Número de passos de tempo: {len(times)}")

    # Plot the transient evolution of Ez at the monitored point
    sa.plot_numerical_time_evolution(times_snap, ez_snap, driver)

    return times, ez_values
    

def plot_fields(total_t, ez_total, incident_t, ez_incident, point=(1, 0)):
    """
    Calcula e plota o coeficiente de reflexão a partir de dados de campo incidentes e totais,
    lidando com possíveis diferenças nos vetores de tempo por meio de interpolação.

    Args:
        total_t (list or np.ndarray): Vetor de tempo para a simulação com PML.
        ez_total (list or np.ndarray): Sinal de campo total (incidente + refletido).
        incident_t (list or np.ndarray): Vetor de tempo para a simulação de referência.
        ez_incident (list or np.ndarray): Sinal de campo incidente da simulação de referência.
    """
    # Garante que todas as entradas sejam arrays NumPy
    total_t = np.array(total_t)
    ez_total = np.array(ez_total)
    incident_t = np.array(incident_t)
    ez_incident = np.array(ez_incident)

    plt.figure(figsize=(12, 6))    
    plt.plot(total_t, ez_total, label='$E_{total}(t)$')
    plt.plot(incident_t, ez_incident, label='$E_{incidente}(t)$', linestyle='--')
    #plt.plot(total_t, ez_reflected, label='$E_{refletida}(t)$', alpha=0.7)
    plt.xlabel("Tempo (s)")
    plt.ylabel("$E_z$ (V/m)")
    plt.title(f"Regime Transitório do Campo Elétrico $E_z$ no Ponto (x={point[0]}, y={point[1]})")
    plt.legend()
    plt.grid(True)


def main() -> None:
    """Função principal para execução do script."""
    clear_terminal() 
           
    # --- Calcular a solução numérica do campo incidente no ponto (x=1, y=0) ---
    mesh_domain, sp = maxwell_at_free_space(PROBLEM, Lx=12.0)
    incident_t, ez_incident = transient_analysis(sp, mesh_domain)

    # --- Calcular a solução numérica do campo total no ponto (x=1, y=0) ---
    mesh_domain, sp = maxwell_at_pml_space(PROBLEM, L=1.0)
    total_t, ez_total = transient_analysis(sp, mesh_domain)

    # --- Calcular e plotar o coeficiente de reflexão ---
    plot_fields(total_t, ez_total, incident_t, ez_incident)
    plt.show()  # Exibe todos os gráficos pendentes


if __name__ == '__main__':
    main()
