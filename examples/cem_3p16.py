#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# pylint: disable=E0611, C0413, C0103, W0401, W0614, W0622
"""
══════════════════════════════════════════════════════════════════════════

ESTE SCRIPT UTILIZA CAMINHOS ABSOLUTOS E NÃO ALTERA O DIRETÓRIO DE TRABALHO

EXECUÇÃO:
cd C:\\git\\PyDG1D
conda activate pyDG1D
python examples\\cem_3p16.py

══════════════════════════════════════════════════════════════════════════

O script está estruturado para garantir que todos os arquivos e pastas usados
sejam acessados por caminhos absolutos, evitando problemas causados por mudanças
do diretório de trabalho (cwd).

PRINCIPAIS CONSIDERAÇÕES SOBRE DIRETÓRIOS:

1️⃣ Diretório raiz do projeto (CWD_ROOT):
    - CWD_ROOT = Path.cwd()
    - Definido automaticamente como o diretório atual quando o script é iniciado.
    - No seu ambiente, o diretório raiz é sempre:
      C:\\Users\\adilt\\OneDrive\\05_GIT\\PYDG1D
    - Todas as pastas e arquivos são referenciados em relação a esse diretório.

2️⃣ Estrutura esperada:
    - PYDG1D/
        ├── examples/
        │   └── jacobi_poly.py
        │   └── vandermonde_matrices.py
        │   └── Cavity1D.py
        │   └── cem_3p16.py
        │   └── hesthaven_e24.py
        │   └── ...
        │   └── ProblemSet1.py
        │   └── LinAdvecEq1D.py
        │   └── LinAdvecEq1D.ipynb
        ├── examplesData/
        │   └── inputs/
        │       ├── LinAdvecEq1D
        │           └── LinAdvecEq1D.mat
        │       ├── jacobi_poly
        │           └── ...
        │       ├── vandermonde_matrices
        │           └── ...
        │   └── outputs/
        │       ├── LinAdvecEq1D
        │           └── LinAdvecEq1D.log
        │       ├── ProblemSet1
        │           └── ProblemSet1.log
        ├── maxwell/
        │   └── dg/
        │       ├── __init__.py
        │       ├── dg1d_tools.py
        │           └── jacobi_polynomial()
        │           └── jacobiGL()

Autor: Adilton Pereira
Data: 26/05/2025
"""

import os
import sys

# Adiciona a raiz do projeto ao PYTHONPATH
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from pathlib import Path
from datetime import datetime
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from maxwell.dg.mesh1d import *
from maxwell.dg.dg1d import *
from maxwell.integrators.LSERK4 import *
from maxwell.utils import *


def resonant_cavity_field(problem_data, x, t):
    """
    Calcula a solução analítica u(x, t) = sin(kx/lambda * (x - a t)).

    Parâmetros
    ----------
    problem_data : dict
        Dados do problema contendo 'kx', 'lmbda', 'advec_speed' e 'tmax'.
    x : ndarray
        Coordenadas físicas onde a solução será avaliada.

    Retorno
    -------
    ndarray
        Valores da solução analítica nos pontos x.
    """
    kx = problem_data['kx']
    lamb = problem_data['lmbda']
    a = problem_data['advec_speed']
    return np.sin(kx / lamb * (x - a * t))


class LinAdvecDriver1D:
    """
    Classe para simulação da equação de advecção linear 1D usando
    o método dos elementos descontínuos (Discontinuous Galerkin - DG).
    A classe implementa um integrador temporal de Runge-Kutta de 5ª ordem
    e calcula o lado direito da equação de advecção.
    """
    def __init__(self, problem_data, sp):
        self.problem_data = problem_data        # Dados do problema
        self.sp = sp                            # Objeto de discretização espacial DG1D
        self.a = problem_data['advec_speed']    # Velocidade de advecção
        CFL = problem_data['cfl']               # Número de Courant-Friedrichs-Lewy
        
        # Condição inicial u(x, 0)
        self.u_h = resonant_cavity_field(problem_data, sp.x, 0.0) 
        self.energy_history = []
        self.time_history = []
        
        # Cálculo do passo de tempo
        xmin = np.min(np.abs(sp.x[0, :] - sp.x[1, :]))
        self.dt = 0.5 * CFL / (2 * np.pi) * xmin

    def AdvecRHS1D(self, u, t, a):
        """
        Computa o lado direito (RHS) da equação de advecção linear 1D
        no contexto do método dos elementos descontínuos (Discontinuous Galerkin).

        Esta função implementa a formulação semi-discreta da equação de advecção,
        incluindo cálculo dos fluxos numéricos nas interfaces dos elementos e
        aplicação das condições de contorno adequadas.

        Parâmetros
        ----------
        u : ndarray of shape (Np, K)
            Solução atual da variável escalar u(x, t), definida nos pontos de colocação
            de cada elemento (Np: pontos por elemento, K: número total de elementos).

        t : float
            Tempo atual da simulação.

        a : float
            Velocidade da advecção. Pode ser positiva ou negativa.

        Retorna
        -------
        rhsu : ndarray of shape (Np, K)
            Termo do lado direito da equação semidiscreta (RHS), pronto para ser usado
            em métodos de integração temporal como Runge-Kutta.

        Notas
        -----
        - Utiliza fluxo upwind com alpha = 1.
        - O fluxo nas interfaces é dado por:
            du = (u⁻ - u⁺) * (a·n - (1 - alpha)|a·n|) / 2
        - As condições de contorno são:
            * Entrada (x = 0): condição senoidal u_in = -sin(a·t)
            * Saída  (x = L): condição de Neumann homogênea (du = 0)

        Etapas
        ------
        1. Achatamento de u e nx no estilo MATLAB (ordem de coluna).
        2. Cálculo do fluxo numérico du nas interfaces.
        3. Substituição de du nas bordas conforme condição de contorno.
        4. Reconstrução de du em forma matricial.
        5. Cálculo final do RHS combinando derivada espacial e termos de contorno.
        """
        sp = self.sp
        if sp.fluxType == "Upwind":
            alpha = 0.0
        elif sp.fluxType == "Centered":
            alpha = 1.0
        else:
            raise ValueError(f"Tipo de fluxo desconhecido: {sp.fluxType}")

        K = sp.mesh.number_of_elements()

        # Flatten estilo MATLAB (coluna a coluna)
        u_flat = u.flatten(order='F')
        nx_flat = sp.nx.flatten(order='F')

        # Condições de contorno Periódicas
        sp.vmap_p[0] = sp.vmap_m[-1]
        sp.vmap_p[-1] = sp.vmap_m[0]

        # Fluxo numérico nas interfaces (du tem shape (2*K,))
        delta_u = u_flat[sp.vmap_m] - u_flat[sp.vmap_p]
        du_flat =  delta_u * (a * nx_flat - (1 - alpha) * np.abs(a * nx_flat)) / 2.0

        # # Condições de contorno
        # uin = -np.sin(a * t)
        # du_flat[sp.map_b[0]] = (u_flat[sp.vmap_b[0]] - uin) * (
        #     a * nx_flat[sp.map_b[0]] -
        #     (1 - alpha) * np.abs(a * nx_flat[sp.map_b[0]])
        # ) / 2.0
        # du_flat[sp.map_b[1]] = 0.0

        # Redimensiona de volta para (2, K)
        du = du_flat.reshape((sp.n_faces, K), order='F')

        # RHS da equação semidiscreta
        rhsu = -a * sp.rx * (sp.diff_matrix @ u) + sp.lift @ (sp.f_scale * du)
        return rhsu

    def run(self, FinalTime):
        """ Executa a simulação e registra energia a cada passo de tempo. """
        time = 0.0
        u = self.u_h.copy()
        resu = np.zeros_like(u)
        Nsteps = int(np.ceil(FinalTime / self.dt))
        self.dt = FinalTime / Nsteps        

        for tstep in range(Nsteps):
            for INTRK in range(5):
                tloc = time + LSERK4.C[INTRK] * self.dt
                rhsu = self.AdvecRHS1D(u, tloc, self.a)
                resu = LSERK4.A[INTRK] * resu + self.dt * rhsu
                u = u + LSERK4.B[INTRK] * resu
            time += self.dt

            # Salvar energia e tempo
            self.energy_history.append(self.ComputeDiscreteEnergy(u))
            self.time_history.append(time)

        self.u_h = u
        return u

    def DriverLog(self, folder_name) -> None:
        """
        Gera um log completo da simulação AdvecDriver1D, salvando em disco e
        exibindo no terminal. Inclui informações sobre discretização espacial,
        inicialização StartUp1D e estado do driver.

        Parâmetros
        ----------
        sp : DG1D
            Objeto de discretização espacial.
        driver : AdvecDriver1D
            Driver da simulação contendo parâmetros e solução numérica.
        output_dir : Path
            Caminho absoluto da pasta de saída onde o log será salvo.
        """
        timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        ua = resonant_cavity_field(self.problem_data, self.sp.x, self.problem_data['tmax'])
        R = np.corrcoef(ua, self.u_h)

        # 1. Informações de discretização espacial
        spatial_info = [
            "\n1. Informações de Discretização Espacial",
            "-" * 75,
            f"Dimensão espacial: {self.sp.mesh.dimension}",
            f"Ordem de interpolação polinomial, N: {self.sp.n_order}",
            f"Número de elementos, K: {self.sp.mesh.number_of_elements()}",
            f"Vértices, VX: {self.sp.mesh.vx}",
            f"Número de vértices: {self.sp.mesh.vx.shape[0]}",
            f"EToV:\n{self.sp.mesh.EToV}"
        ]

        # 2. StartUp1D
        startup_info = [
            "\n2. Inicialização - StartUp1D",
            "-" * 75,
            f"Pontos de colocação (LGL): {self.sp.r}",
            f"\nCoordenadas nodais x:\n{self.sp.x}",
            f"\nJacobiano:\n{self.sp.jacobian}",
            f"\nJacobiano inverso rx:\n{self.sp.rx}",
            #f"\nMatriz de conectividade EToE:\n{self.sp.etoe}",
            #f"\nMatriz de conectividade EToF:\n{self.sp.etof}",
            "\nMapas de conectividade:",
            f"  vmapM: {self.sp.vmap_m}",
            f"  vmapP: {self.sp.vmap_p}",
            f"  vmapB: {self.sp.vmap_b}",
            f"   mapB: {self.sp.map_b}",
            f"\nNormais nas interfaces nx:\n{self.sp.nx}"
        ]

        # 3. Informações do driver
        driver_info = [
            "\n3. Informações do Driver",
            "-" * 75,
            f"Tipo de fluxo: {self.sp.fluxType}",
            f"Velocidade de advecção, a: {self.a}",
            f"Número de Courant-Friedrichs-Lewy, CFL: {self.problem_data['cfl']}",
            f"Passo de tempo, dt: {self.dt:.4e}",
            f"Tempo final, FinalTime: {self.problem_data['tmax']:.4e}",
            f"Condição de contorno: {self.sp.mesh.boundary_label}",

            format_matrix(self.u_h, title=f"DG Solution u_h @ t = {self.problem_data['tmax']:.2f}"),
            format_matrix(ua, title=f"Solução Analítica u_a @ t = {self.problem_data['tmax']:.2f}"),            
            f"\n🌐 Correlação entre solução analítica e numérica: {R[0, 1]:.2f}"
        ]

        header = [
            "=" * 80,
            "AdvecDriver1D - Log de Execução",
            f"Data e Hora: {timestamp}",
            "=" * 80
        ]

        footer = ["=" * 80]

        # Combina todas as seções
        log_lines = header + spatial_info + startup_info + driver_info + footer

        # Criação da pasta de saída, se necessário
        # INPUTS = (Path(__file__).parent.parent / 'examplesData' / 'inputs' / PROBLEM['name']).resolve()
        OUTPUTS = (Path(__file__).parent.parent / 'examplesData' / 'outputs' / folder_name).resolve()
        OUTPUTS.mkdir(parents=True, exist_ok=True)
        log_path = OUTPUTS / f"{self.problem_data['name']}.log"

        # Salvando no arquivo
        with open(log_path, 'w', encoding='utf-8') as f:
            f.write('\n'.join(log_lines))

        # Impressão no terminal
        print("\n".join(log_lines))
        print("\n📄 LOG GERADO:", log_path)

    def ComputeDiscreteEnergy(self, u):
        """
        Calcula a energia discreta associada à solução numérica u.

        Parâmetros
        ----------
        u : ndarray of shape (Np, K)
            Solução DG em todos os elementos.

        Retorno
        -------
        float
            Valor escalar da energia no tempo atual.
        """
        M = self.sp.mass
        J = self.sp.jacobian
        energy = 0.0
        for k in range(self.sp.mesh.number_of_elements()):
            uk = u[:, k][:, np.newaxis]
            Jk = np.diag(J[:, k])
            energy += (uk.T @ Jk @ M @ uk)[0, 0]
        return energy

    def ComputeL2Error(self, ua):
        """
        Calcula o erro L2 entre a solução numérica u_h e a solução analítica ua.
        Parâmetros
        ----------
        sp : DG1D
            Objeto de discretização espacial.
        ua : ndarray
            Solução analítica avaliada nos pontos nodais.
        Retorno
        -------
        float
            Erro L2 global entre a solução numérica e a analítica.
        """
        
        err = ua - self.u_h
        M = self.sp.mass
        J = self.sp.jacobian
        K = self.sp.mesh.number_of_elements()
        errL2_local = np.zeros(K)

        for k in range(K):
            Jk = np.diag(J[:, k])
            ek = err[:, k][:, np.newaxis]
            errL2_local[k] = (ek.T @ Jk @ M @ ek)[0, 0]

        return np.sqrt(np.sum(errL2_local))


def plot_solution(problem, sp, u_h) -> None:
    """
    Retorna a solução analítica u(x, t) = sin(x - a t) e
    a solução numérica u_h(x, t) em um gráfico.

    Parâmetros
    ----------
    sp : DG1D
        Objeto de discretização espacial.
    a : float
        Velocidade de advecção.
    FinalTime : float
        Tempo final da simulação.

    Retorna
    -------
    np.ndarray
        Vetor com a solução analítica nos pontos nodais.
    """
    K = sp.mesh.number_of_elements()
    xmin, xmax = sp.mesh.xmin, sp.mesh.xmax
    x_scale = problem['L']
    xa = np.linspace(xmin, xmax, 1000)
    ua = resonant_cavity_field(problem, xa, problem['tmax'])
    _, axs = plt.subplots(1, 1, figsize=(10, 6), sharey=True)

    # === Subplot 1: Solução global ===
    for k in range(K):
        axs.plot(
            sp.x[:, k] / x_scale, u_h[:, k],
            marker='o',
            markerfacecolor='none',
            markersize=3,
            linestyle='-',
            linewidth=1.5,
            label = fr'$D^{{{k+1}}}$'
        )

    # Sobreposição da solução analítica
    axs.plot(xa / x_scale, ua, 'k--', linewidth=1, label='Solução analítica')
    axs.set_title(f"1D linear advection equation with {sp.fluxType} Flux"
        f"\n a = {problem['advec_speed']:.2f}, N = {problem['n_order']}, "
        f"K = {K}, T = {problem['tmax']:.0f}")
    axs.set_xlabel(fr"$x/{x_scale:.2f}$")
    axs.set_ylabel(r"$u(x)$")
    axs.set_xlim(xmin, xmax / x_scale)
    axs.set_ylim(-1.5, 1.5)
    axs.set_xticks(np.linspace(xmin, xmax / x_scale, K+1))
    axs.legend(
        loc='center left',
        bbox_to_anchor=(1.02, 0.50),
        ncol=2,
        fontsize=9,
        frameon=False
    )

    plt.grid(False)
    plt.tight_layout()
    plt.show()


def plot_energy_history(driver) -> None:
    """
    Plota a evolução temporal da energia da solução DG.
    """
    plt.figure(figsize=(8, 4))
    plt.plot(driver.time_history, driver.energy_history, 'b-', linewidth=2)
    plt.title("Evolução temporal da energia discreta")
    plt.xlabel("Tempo")
    plt.ylabel("Energia")
    plt.grid(True)
    plt.tight_layout()
    plt.show()


def single_test_solution(problem) -> None:
    """
    Testa a solução numérica comparando com a solução analítica.

    Parâmetros
    ----------
    problem : dict
        Dicionário com os parâmetros do problema.
    sp : DG1D
        Objeto de discretização espacial.
    uh : ndarray
        Solução numérica obtida pelo driver.

    Retorna
    -------
    None
    """
    
    # Generate Spatial Discretization
    sp = DG1D(
        n_order=problem['n_order'],
        mesh=Mesh1D(
            xmin=0.0,
            xmax=problem['L'],
            k_elem=problem['k_elem'],
            boundary_label="Periodic"),
        fluxType=problem['flux_type'])
    
    # Initialize the solver
    driver = LinAdvecDriver1D(problem, sp)
    
    # Set the time integrator
    driver.run(FinalTime=problem['tmax'])

    # Logging
    driver.DriverLog(folder_name='cem_3p1')

    # # Compute the correlation coefficient
    # u_expected = resonant_cavity_field(problem, sp.x, problem['tmax'])
    # print(f"\n🌐 Solução esperada (analítica) em t = {problem['tmax']}:")
    # print(f"\n{u_expected}")
    # R = np.corrcoef(u_expected, driver.u_h)
    # print(f"\n🌐 Correlação entre solução esperada e solução numérica: {R[0,1]:.2f}")

    # Compute L2 error
    # ua = analytical_solution(problem, sp.x)
    ua = resonant_cavity_field(problem, sp.x, problem['tmax'])
    L2_error = driver.ComputeL2Error(ua)
    print(f"\n🌐 Erro global na norma L2: {L2_error:.1e}")

    # Plotting the solution
    plot_solution(problem, sp, driver.u_h)

    # Plotting discrete energy
    # driver.plot_energy_history()


def run_convergence(problem):
    """ 
    Executa um estudo de convergência para a equação de advecção linear 1D.
    Parâmetros
    ----------
    N : int
        Ordem de interpolação polinomial para o método DG1D.
    
    Retorno
    -------
    pd.DataFrame
    DataFrame contendo os erros L2 e taxas de convergência para diferentes
    números de elementos (K).
    """
    errors = []
    rates = [None]
    k_list = [2, 4, 8, 16, 32, 64]
    for k in k_list:
        # Generate Spatial Discretization
        sp = DG1D(
            n_order=problem['n_order'],
            mesh=Mesh1D(
                xmin=0.0,
                xmax=problem['L'],
                k_elem=k,
                boundary_label="Periodic"),
            fluxType=problem['flux_type'])
        
        # Initialize the solver
        driver = LinAdvecDriver1D(problem, sp)
        
        # Compute L2 error
        driver.run(FinalTime=problem['tmax'])
        ua = resonant_cavity_field(problem, sp.x, problem['tmax'])
        errors.append(driver.ComputeL2Error(ua))

    for i in range(1, len(errors)):
        if errors[i] == 0 or errors[i-1] == 0:
            rates.append(None)
        else:
            rates.append(np.log(errors[i-1] / errors[i]) / np.log(2))

    df = pd.DataFrame({
        "K": k_list,
        "Erro L2": errors,
        "Taxa de Convergência": rates
    })

    pd.set_option('display.precision', 2)
    pd.set_option('expand_frame_repr', False)
    print(f"\nTabela de Erros e Taxas de Convergência (Ordem N = {problem['n_order']})")
    print(df.to_string(index=False))
    return df


def compute_energy_evolution(problem: dict):
    """
    Executa simulações para diferentes valores de K e plota a evolução
    da energia discreta. Retorna também as soluções numéricas e malhas
    para uso posterior.

    Parâmetros
    ----------
    problem : dict
        Dicionário com os parâmetros do problema.

    Retorno
    -------
    results : list of dict
        Lista contendo, para cada K:
        {
            'K': número de elementos,
            'x': coordenadas nodais (Np x K),
            'u_h': solução final DG (Np x K),
            'sp': objeto DG1D (opcional, se quiser usar em outros contextos),
            'driver': objeto LinAdvecDriver1D
        }
    """
    k_list = [2, 4, 8, 16, 32, 64]
    k_list = [2, 4, 8]
    plt.figure(figsize=(10, 6))
    results = []

    for K in k_list:
        problem['k_elem'] = K

        sp = DG1D(
            n_order=problem['n_order'],
            mesh=Mesh1D(
                xmin=0.0,
                xmax=problem['L'],
                k_elem=K,
                boundary_label="Periodic"),
            fluxType=problem['flux_type'])

        driver = LinAdvecDriver1D(problem, sp)
        driver.run(FinalTime=problem['tmax'])

        # Plot da energia
        plt.plot(driver.time_history, driver.energy_history, label=fr'$K = {K}$')

        # Armazena resultados para uso posterior
        results.append({
            'K': K,
            'x': sp.x.copy(),
            'u_h': driver.u_h.copy(),
            'sp': sp,
            'driver': driver
        })

    plt.title("Evolução temporal da energia discreta para diferentes K")
    plt.xlabel("Tempo")
    plt.ylabel("Energia")
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.show()

    return results


def plot_multi_solution(problem: dict, results: list) -> None:
    """
    Gera subplots da solução numérica final para diferentes K usando
    os dados já computados e armazenados em `compute_energy_evolution`.

    Parâmetros
    ----------
    results : list of dict
        Lista de dicionários com as chaves 'K', 'x', 'u_h' retornadas pela
        função compute_energy_evolution.

    problem : dict
        Parâmetros do problema, usado para gerar a solução analítica.
    """
    nrows, ncols = 2, 3
    fig, axs = plt.subplots(nrows, ncols, figsize=(15, 8), sharey=True)
    axs = axs.flatten()

    for idx, data in enumerate(results):
        K = data['K']
        x = data['x']
        u_h = data['u_h']
        xa = np.linspace(0, problem['L'], 1000)
        ua = resonant_cavity_field(problem, xa, problem['tmax'])

        ax = axs[idx]
        for k in range(K):
            ax.plot(
                x[:, k],
                u_h[:, k],
                marker='o',
                markerfacecolor='none',
                markersize=3,
                linestyle='-',
                linewidth=1.0
            )
        ax.plot(xa, ua, 'k--', linewidth=1.2, label="Analítica")
        ax.set_title(fr'$K = {K}$')
        ax.set_xlim([0, problem['L']])
        ax.set_ylim([-1.5, 1.5])
        ax.grid(False)
        ax.set_xlabel('$x$')
        if idx % ncols == 0:
            ax.set_ylabel('$u(x)$')

    fig.suptitle("Soluções numéricas $u_h(x)$ e analítica para diferentes K", fontsize=14)
    plt.tight_layout(rect=[0, 0, 1, 0.95])
    plt.show()


def main() -> None:
    """Função principal para execução do script."""
    clear_terminal()

    # PROBLEMA: Definições do problema de advecção linear 1D
    PROBLEM = {'name': 'cem_3p16',
        'description': 'Problem 1 - Computational Electromagnetics List 3',
        'advec_speed': 1,           # Velocidade de advecção
        'flux_type': 'Upwind',      # 'Upwind' or 'Centered'
        'cfl': 0.75,                # Número de Courant-Friedrichs-Lewy
        'tmax': 1E1,                # Tempo final da simulação
        'kx': 2*np.pi,              # Número de onda
        'lmbda': 1,                 # Comprimento de onda
        'L': 1,                     # Comprimento do domínio
        'n_order': 2,               # Ordem de interpolação polinomial
        'k_elem': 8                 # Número de elementos na malha
    }

    # Test 0 - The solution
    single_test_solution(PROBLEM)

    # # Test 1 - Run convergence study
    # run_convergence(PROBLEM)

    # # Test 2 - Compute energy evolution
    # results = compute_energy_evolution(PROBLEM)
    # plot_multi_solution(PROBLEM, results)


if __name__ == '__main__':
    main()
