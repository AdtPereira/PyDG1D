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


PROBLEM = {
    'DIM': 2,                       # Dimensão do problema (2D)
    'name': 'single_conductor',
    'folder': 'single_conductor',
    'description': 'Análise espectral de um pulso gaussiano propagante em 2D',
    'num_snapshots': 9,             # Número de snapshots a serem coletados
    'dg': {
        'n_order': 3,               # Ordem de interpolação polinomial
        'flux_type': 'Upwind',      # 'Upwind' or 'Centered'
        'cfl': 1.0,                 # Número de Courant-Friedrichs-Lewy
        'bc': "PEC",                # Condição de contorno: 'PEC', 'SMA'  or 'Periodic
        't_final': 3.0,             # Tempo final da simulação
    },
    'domain': {
        'type': 'rectangle',    # Tipo de domínio: 'rectangle' ou 'circle'
        'Lx': 2.0,              # Largura do domínio físico na direção x
        'Ly': 2.0,              # Largura do domínio físico na direção y
        'rc': 0.5               # Raio do círculo central (condutor)
    },
    'source': {
        'type': 'gaussian', # Tipo de fonte: 'gaussian' ou 'derivative'
        'A0': 1.0,          # Amplitude do pulso gaussiano
        'sigma_r': 0.2,     # Desvio padrão espacial do pulso gaussiano
        'sigma_t': 2e-7,    # Desvio padrão temporal do pulso gaussiano
        't0': 1e-6          # Tempo de pico do pulso gaussiano
    },
    'pml': {
        'type': 'Abarbanel and Gottlieb', 
        'x0': 1.0,          # Fronteira do domínio físico
        'y0': 1.0,          # Fronteira do domínio físico
        'L': 1.0,           # Largura da camada da PML
        'pml_order': 2,     # Ordem polinomial da PML
        'R': 1E-4           # Coeficiente de reflexão na interface da PML
    }
}


class SpectralAnalyzer:
    def __init__(self, problem, mesh):
        self.N = problem['dg']['n_order']               # Ordem de interpolação polinomial
        self.L = 2*problem['domain']['Lx']              # Largura do domínio físico
        self.t_final = problem['dg']['t_final']         # Tempo final da simulação
        self.num_snapshots = problem['num_snapshots']   # Número de snapshots a serem coletados

        self.t_list = np.linspace(0, self.t_final, self.num_snapshots)
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


    def plot_source_spectrum(self):
        """
        Avalia e plota o espectro de potência em dB, otimizado para a faixa de MHz.

        Esta função usa parâmetros em microssegundos (µs) e exibe os eixos dos gráficos
        em microssegundos e Megahertz (MHz) para uma análise clara de sinais de
        baixa frequência.

        Parâmetros:
        -----------
        sigma_t_us : float, opcional
            Desvio padrão temporal (largura) do pulso em MICROSSEGUNDOS (µs).
            Default: 0.2.
        simulation_time_us : float, opcional
            Duração total da simulação no domínio do tempo, em MICROSSEGUNDOS (µs).
            Default: 6.0.
        num_samples : int, opcional
            Número de amostras para a série temporal. Default: 2**14.
        """
        sigma_t = PROBLEM['source']['sigma_t']  # Desvio padrão temporal do pulso gaussiano
        t_max = PROBLEM['dg']['t_final']  # Tempo final da simulação
        num_samples = 2**14         # Número de amostras para a série temporal

        # 2. Configurar o vetor de tempo para a análise
        time_vector = np.linspace(0, t_max, num_samples, endpoint=False)
        dt = time_vector[1] - time_vector[0]

        # 3. Gerar o pulso gaussiano no domínio do tempo
        pulse_in_time = self.gaussian_pulse(x=0, y=0, t=time_vector)

        # 4. Calcular a FFT e as frequências correspondentes em Hz
        pulse_in_freq = np.fft.fft(pulse_in_time)
        frequencies_hz = np.fft.fftfreq(num_samples, d=dt)

        # 5. Focar apenas nas frequências positivas
        positive_freq_mask = frequencies_hz >= 0
        freqs_positive_hz = frequencies_hz[positive_freq_mask]
        magnitude_spectrum = np.abs(pulse_in_freq[positive_freq_mask])

        # 6. Converter a magnitude para potência em decibéis (dB)
        epsilon = 1e-20
        max_magnitude = np.max(magnitude_spectrum)
        if max_magnitude > epsilon:
            power_db = 20 * np.log10(magnitude_spectrum / max_magnitude + epsilon)
        else:
            power_db = np.full_like(magnitude_spectrum, -np.inf)

        # --- Plotagem dos Resultados ---
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6))
        fig.suptitle(fr'Gaussian $E_z(z,y,t)$ Source ($\sigma_t={1e6*sigma_t:.1f}$ µs, $t_{{max}}$={1e6*t_max} µs)', fontsize=14)

        # Gráfico do Domínio do Tempo (eixo em µs)
        ax1.plot(time_vector * 1e6, pulse_in_time, color='red')
        ax1.set_title("Pulso Gaussiano no Domínio do Tempo")
        ax1.set_xlabel("Tempo (µs)")
        ax1.set_ylabel("Amplitude")
        ax1.grid(True, linestyle=':')
        ax1.set_xlim(0, 1e6*t_max)

        # Gráfico do Espectro de Potência (eixo em MHz)
        ax2.plot(freqs_positive_hz / 1e6, power_db, color='red', marker='.', markersize=4, linestyle='-')
        ax2.set_title("Espectro de Potência")
        ax2.set_xlabel("Frequência (MHz)")
        ax2.set_ylabel("Potência (dB)")
        ax2.grid(True, linestyle=':')

        # Ajusta o limite do eixo de frequência para a faixa de MHz
        freq_cutoff_mhz = (1.5 / sigma_t) / 1e6
        ax2.set_xlim(0, freq_cutoff_mhz)
        ax2.set_ylim(-70, 5)
        plt.tight_layout(rect=[0, 0.03, 1, 0.95])


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
            plt.scatter(xk, yk, marker=mk, facecolors='none', color=ck, s=3, label=f'e{k}')

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


    def plot_ez_transient_analysis(self, times, ez_values):
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


    def plot_field_evolution_on_x_axis(self, sa, sp, x0, N_grid=500):
        """
        Monitora e plota o campo Ez ao longo do eixo x (y=0) em vários instantes de tempo,
        com cada instante exibido em um subplot separado.

        Parâmetros:
        -----------
        SpectralAnalyzer : object
            Instância da classe que contém os métodos e parâmetros do problema.
        sp : Maxwell2D
            Objeto DG contendo a malha e as coordenadas dos pontos.
        N_grid : int
            Número de pontos ao longo do eixo x para a interpolação.
        """
        print("\n🚀 Iniciando simulação para monitorar Ez no eixo x (em subplots)...")

        # --- Fase de Simulação ---
        driver = MaxwellDriver(sp, CFL=PROBLEM['dg']['cfl'])
        driver['Ez'][:] = sa.gaussian_pulse(sp.x, sp.y, 0)
        # driver['Ez'][:] = SpectralAnalyzer.gaussian_pulse_derivative(sp.x, sp.y, 0)

        print(f"\n🌐 Inicializando o driver Maxwell com CFL = {PROBLEM['dg']['cfl']:.2f} e dt = {1E3*driver.dt:.1f} ms")

        ez_snapshots = []
        for t in self.t_list:
            driver.run_until(t)
            ez_snapshots.append(driver['Ez'].copy())
            print(f"  -> Snapshot capturado em t = {1E3*t:.0f} ms")

        # --- Fase de Extração e Pré-cálculo ---
        x_nodes = driver.sp.x.ravel(order='F')
        y_nodes = driver.sp.y.ravel(order='F')
        x_min, x_max = x_nodes.min(), x_nodes.max()
        x_line = np.linspace(x_min, x_max, N_grid)
        points_on_x_axis = np.column_stack((x_line, np.zeros_like(x_line)))

        # Interpola todos os dados primeiro para encontrar os limites globais do eixo y.
        # Isso garante que todos os subplots tenham a mesma escala vertical.
        all_ez_on_axis = []
        for Ez_snapshot in ez_snapshots:
            u_snapshot = Ez_snapshot.ravel(order='F')
            ez_on_axis = griddata(points=(x_nodes, y_nodes),
                                  values=u_snapshot,
                                  xi=points_on_x_axis,
                                  method='cubic')
            all_ez_on_axis.append(ez_on_axis)
        
        # Define os limites verticais (y-axis)
        vmin = np.nanmin([np.nanmin(data) for data in all_ez_on_axis])
        vmax = np.nanmax([np.nanmax(data) for data in all_ez_on_axis])
        padding = (vmax - vmin) * 0.05 # Adiciona um pequeno preenchimento
        y_limits = (vmin - padding, vmax + padding)


        # --- Fase de Plotagem ---
        ncols = int(np.ceil(np.sqrt(self.num_snapshots)))
        nrows = int(np.ceil(self.num_snapshots / ncols))

        fig, axes = plt.subplots(nrows, ncols, figsize=(5 * ncols, 3.5 * nrows), 
                                 constrained_layout=True, sharex=True, sharey=True)
        axes = axes.ravel()

        fig.suptitle(r'Evolução do Campo $E_z$ ao longo do eixo x ($y=0$)', fontsize=16)

        for i, (ez_data, t) in enumerate(zip(all_ez_on_axis, self.t_list)):
            ax = axes[i]

            # >>> NOVIDADE: Colorir as regiões da PML <<<
            # Colore a PML da direita (x > interface)
            ax.axvspan(x0, x_max, color='gray', alpha=0.3, zorder=0)
            # Colore a PML da esquerda (x < -interface)
            ax.axvspan(x_min, -x0, color='gray', alpha=0.3, zorder=0)

            ax.plot(x_line, ez_data, color='dodgerblue')
            ax.set_title(f"t = {1E3*t:.0f} ms")
            ax.grid(True, linestyle='--', alpha=0.6)
            ax.set_ylim(y_limits) # Aplica os mesmos limites a todos os subplots

        # Configura rótulos comuns para a figura inteira
        fig.supxlabel('Posição x (m)')
        fig.supylabel(r'$E_z$ (V/m)')

        # Desativa eixos extras que não estão sendo usados
        for j in range(i + 1, len(axes)):
            axes[j].axis('off')
        
        print("\n✅ Gráfico com subplots gerado com sucesso.")


def plotar_condutividade_por_elemento(sp, ax=None, **kwargs):
    """
    Plota um mapa de cores do perfil de condutividade definido por elemento (sigma_c).

    A função mapeia o valor de condutividade de cada elemento para todos os nós
    de visualização dentro dele, criando um mapa de cores preciso da propriedade
    constante por elemento.

    Parâmetros:
    -----------
    sp : object
        O objeto de espacialização (ex: Maxwell2D) que contém a malha e os dados.
    ax : matplotlib.axes.Axes, opcional
        O eixo onde o gráfico será desenhado. Se None, cria uma nova figura.
    **kwargs :
        Argumentos extras para a função `tricontourf` (ex: cmap, levels).
    """
    # Verifica se há um eixo, se não, cria um novo
    show_plot = False
    if ax is None:
        fig, ax = plt.subplots(figsize=(8, 7))
        show_plot = True

    # Nível de refinamento para a visualização
    Nout = sp.n_order

    # O vetor sp.sigma (sigma_c) tem formato (K,). Precisamos expandi-lo para o
    # formato de campo (Np, K) para usar a rotina de plotagem.
    Np = sp.number_of_nodes_per_element()
    # np.tile repete o vetor de condutividade para cada nó do elemento
    sigma_field = np.tile(sp.sigma, (Np, 1))

    # --- Lógica de visualização baseada em `plot_field` ---
    # (Esta é a maneira correta de visualizar campos em malhas DG)

    # Constrói grade igualmente espaçada no triângulo de referência
    # (assumindo que as funções auxiliares estão disponíveis)
    try:
        r, s = xy_to_rs(*set_nodes_in_equilateral_triangle(sp.n_order))
        rout, sout = xy_to_rs(*set_nodes_in_equilateral_triangle(Nout))
        interp = vandermonde(sp.n_order, rout, sout) @ np.linalg.inv(vandermonde(sp.n_order, r, s))
    except (ImportError, NameError):
        # Fallback se as funções não estiverem no escopo - pode ser necessário adaptar
        print("Aviso: Funções de interpolação não encontradas. Usando um fallback simples.")
        # Este fallback é apenas uma simplificação e pode não ser preciso
        rout, sout = np.linspace(-1, 1, 10), np.linspace(-1, 1, 10) # Placeholder
        interp = np.eye(Np) # Placeholder
        Nout = sp.n_order # Placeholder

    # Interpola coordenadas e o campo sigma para os nós de plotagem
    xout = interp @ sp.x
    yout = interp @ sp.y
    uout = interp @ sigma_field  # Usa o campo de sigma que criamos

    # Constrói a triangulação para todos os nós de plotagem
    # Usando Delaunay para robustez
    from scipy.spatial import Delaunay
    tri = Delaunay(np.vstack((rout, sout)).T).simplices
    TRI = np.array([], dtype=int).reshape(0, 3)
    Npout = len(rout)
    for k in range(sp.mesh.number_of_elements()):
        TRI = np.vstack((TRI, tri + k * Npout))

    # Define opções padrão de plotagem
    if 'levels' not in kwargs:
        # Se os valores são discretos, `levels` ajuda a definir as fronteiras
        unique_sigmas = np.unique(uout)
        if len(unique_sigmas) < 20: # Bom para poucos valores de sigma
             kwargs['levels'] = len(unique_sigmas)
        else:
             kwargs['levels'] = 20 # Padrão para muitos valores
    if 'cmap' not in kwargs:
        kwargs['cmap'] = 'turbo' # Um bom mapa de cores para distinguir regiões

    # Renderiza o mapa de cores
    im = ax.tricontourf(
        xout.ravel('F'),
        yout.ravel('F'),
        uout.ravel('F'),
        triangles=TRI,
        **kwargs
    )

    # Adiciona rótulos e formatação
    fig = ax.get_figure()
    fig.colorbar(im, ax=ax, label=r"Condutividade $\sigma_c$ (S/m)")
    ax.set_title(r'Mapa de Condutividade por Elemento ($\sigma_c$)')
    ax.set_xlabel("Coordenada x (m)")
    ax.set_ylabel("Coordenada y (m)")
    ax.axis('equal')
    ax.grid(True, linestyle='--', alpha=0.5)

    if show_plot:
        plt.show()


def aplicar_condutividade(x_nodes, y_nodes, rc, sigma_condutor):
    """
    Aplica uma condutividade física aos elementos cujo baricentro está dentro
    de um condutor circular central.

    Esta função é projetada para métodos de elementos, como Galerkin Discontínuo (DG),
    onde as propriedades são constantes por elemento.

    Parâmetros:
    ----------
    x_nodes : np.ndarray
        Array 2D de formato (Np, K) com as coordenadas x de todos os nós,
        onde Np é o número de nós por elemento e K é o número de elementos.
        (Corresponde a `self.x` na classe Maxwell2D).
    y_nodes : np.ndarray
        Array 2D de formato (Np, K) com as coordenadas y de todos os nós.
        (Corresponde a `self.y` na classe Maxwell2D).
    rc : float
        O raio do condutor circular central.
    sigma_condutor : float
        O valor da condutividade elétrica (σ) a ser aplicado dentro do condutor.

    Retorna:
    -------
    np.ndarray
        Um array 1D de formato (K,), onde K é o número de elementos. Cada entrada
        contém o valor da condutividade para o elemento correspondente.
    """
    # Garante que a entrada tenha o formato esperado (Np, K)
    if x_nodes.ndim != 2 or y_nodes.ndim != 2:
        raise ValueError("As entradas x_nodes e y_nodes devem ser arrays 2D de formato (Nós, Elementos).")

    # Calcula o baricentro (centro geométrico) para cada elemento.
    # A função np.mean com axis=0 calcula a média ao longo do eixo dos nós (Np),
    # resultando em um array de formato (K,).
    x_baricentros = np.mean(x_nodes, axis=0)
    y_baricentros = np.mean(y_nodes, axis=0)

    # O número de elementos é o tamanho do array de baricentros.
    num_elementos = len(x_baricentros)

    # Inicializa o array de condutividade com zeros.
    sigma = np.zeros(num_elementos)

    # Calcula o quadrado da distância de cada baricentro até a origem (0, 0).
    distancia_quadrada = x_baricentros**2 + y_baricentros**2

    # Cria uma máscara booleana para os elementos cujo baricentro está dentro do raio.
    mascara_condutor = distancia_quadrada < rc**2

    # Aplica o valor da condutividade a esses elementos.
    sigma[mascara_condutor] = sigma_condutor

    return sigma


def maxwell_without_conductor(problem, L, h):
    """
    Monitors and plots the electric field Ez at the PML interface (x=1, y=0) over time.
    """
    print(f"\n🌐 Inicializando o analisador espectral para {problem['name']} truncado com PML ...")
    
    pml_project = {
        'x0': 1.0,          # Fronteira direita da PML
        'L': L,             # Largura da camada da PML
        'pml_order': 2,     # Ordem polinomial da PML
        'R': 1E-4,          # Coeficiente de reflexão na interface da PML
    }

    mesh_data = mesh_single_conductor_domain(problem, h, view_mesh=False)
    mesh_domain = Mesh2D(
        vx=mesh_data['VX'],
        vy=mesh_data['VY'],
        EToV=mesh_data['EToV'],
        boundary_label=problem['dg']['bc'])

    sp = Maxwell2D(n_order=problem['dg']['n_order'], mesh=mesh_domain, fluxType=problem['dg']['flux_type'], sigma=None, pml_design=pml_project)
    return mesh_domain, sp


def maxwell_with_conductor(problem, L, h):
    """
    Monitors and plots the electric field Ez at the PML interface (x=1, y=0) over time.
    """
    print(f"\n🌐 Inicializando o analisador espectral para {problem['name']} truncado com PML ...")
    
    pml_project = {
        'x0': 1.0,          # Fronteira direita da PML
        'L': L,             # Largura da camada da PML
        'pml_order': 2,     # Ordem polinomial da PML
        'R': 1E-4,          # Coeficiente de reflexão na interface da PML
    }

    mesh_data = mesh_single_conductor_domain(problem, h, view_mesh=False)
    mesh_domain = Mesh2D(
        vx=mesh_data['VX'],
        vy=mesh_data['VY'],
        EToV=mesh_data['EToV'],
        boundary_label=problem['dg']['bc'])

    # Aplicar condutividade do condutor circular central
    x_nodes, y_nodes = nodes_coordinates(problem['dg']['n_order'], mesh_domain)
    sigma_core = aplicar_condutividade(x_nodes, y_nodes, problem['domain']['rc'], sigma_condutor=1E1)

    sp = Maxwell2D(n_order=problem['dg']['n_order'], mesh=mesh_domain, fluxType=problem['dg']['flux_type'], sigma=sigma_core, pml_design=pml_project)

    # 2. Plotar o mapa de condutividade do condutor (sigma_c)
    plotar_condutividade_por_elemento(sp)


    return mesh_domain, sp


def transient_analysis(sp, mesh_domain, point=(1, 0), n_samples=500, n_snaps=9):
    """
    Monitors and plots the electric field Ez at the PML interface (x=1, y=0) over time.
    """
    print(f"\n🌐 Inicializando o analisador espectral para {PROBLEM['name']} ...")
    sa = SpectralAnalyzer(PROBLEM, mesh_domain)
    sa.plot_collocation_points(show_ids=False, sigma_data=False)
    sa.plot_field_evolution_on_x_axis(sa, sp, x0=1.0, N_grid=500)

    # 1. Inicializar o campo elétrico Ez 
    driver = MaxwellDriver(sp, CFL=PROBLEM['dg']['cfl'])
    driver['Ez'][:] = sa.gaussian_pulse(sp.x, sp.y, 0)

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
           
    # --- Calcular a solução numérica do campo incidente no ponto espaço livre ---
    mesh_domain, sp = maxwell_without_conductor(PROBLEM, L=1.0, h=0.2)
    incident_t, ez_incident = transient_analysis(sp, mesh_domain)

    # --- Calcular a solução numérica do campo total no domínio físico ---
    mesh_domain, sp = maxwell_with_conductor(PROBLEM, L=1.0, h=0.2)
    total_t, ez_total = transient_analysis(sp, mesh_domain)

    # --- Plotar o campo total e o campo incidente ---
    plot_fields(total_t, ez_total, incident_t, ez_incident, point=(1, 0))


if __name__ == '__main__':
    main()
    plt.show()
