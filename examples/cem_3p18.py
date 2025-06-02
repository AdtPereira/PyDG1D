#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# pylint: disable=E0611, C0413, C0103, W0401, W0614, W0622
"""
══════════════════════════════════════════════════════════════════════════

ESTE SCRIPT UTILIZA CAMINHOS ABSOLUTOS E NÃO ALTERA O DIRETÓRIO DE TRABALHO

EXECUÇÃO:
cd C:\\git\\PyDG1D
conda activate pyDG1D
python examples\\cem_3p18.py

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


class LinAdvecDriver1D:
    """
    Classe para simulação da equação de advecção linear 1D usando
    o método dos elementos descontínuos (Discontinuous Galerkin - DG).
    A classe implementa um integrador temporal de Runge-Kutta de 5ª ordem
    e calcula o lado direito da equação de advecção.
    """
    def __init__(self, problem_data, sp, u0):
        self.problem_data = problem_data        # Dados do problema
        self.sp = sp                            # Objeto de discretização espacial DG1D
        self.a = problem_data['advec_speed']    # Velocidade de advecção
        self.u_h = u0                           # Condição inicial u(x, 0)
        CFL = problem_data['cfl']               # Número de Courant-Friedrichs-Lewy
        
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
        kx = self.problem_data['kx']

        # Flatten estilo MATLAB (coluna a coluna)
        u_flat = u.flatten(order='F')
        nx_flat = sp.nx.flatten(order='F')

        # # Condições de contorno Periódicas
        # sp.vmap_p[0] = sp.vmap_m[-1]
        # sp.vmap_p[-1] = sp.vmap_m[0]

        # Fluxo numérico nas interfaces (du tem shape (2*K,))
        delta_u = u_flat[sp.vmap_m] - u_flat[sp.vmap_p]
        du_flat =  delta_u * (a * nx_flat - (1 - alpha) * np.abs(a * nx_flat)) / 2.0

        # Condições de contorno
        uin = -np.sin((kx * a) * t)
        du_flat[sp.map_b[0]] = (u_flat[sp.vmap_b[0]] - uin) * (
            a * nx_flat[sp.map_b[0]] -
            (1 - alpha) * np.abs(a * nx_flat[sp.map_b[0]])
        ) / 2.0
        du_flat[sp.map_b[1]] = 0.0

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
            f"Condição de contorno: {self.sp.mesh.boundary_label}",
            # f"\nSolução numérica u_h:\n{np.array2string(driver.u_h, precision=4)}"
            format_matrix(self.u_h, title="DG Solution u_h")
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
        print("\n📄 LOG GERADO:", log_path)
        print("\n".join(log_lines))

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


def analytical_solution(problem_data, x):
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
    a = problem_data['advec_speed']
    t = problem_data['tmax']
    return -np.sin(kx * a * (t - x/a))


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
    x_ana = np.linspace(xmin, xmax, 1000)
    u_ana = analytical_solution(problem, x_ana)
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
    axs.plot(x_ana / x_scale, u_ana, 'k--', linewidth=1, label='Solução analítica')
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
            k_elem=problem['k_elem']),
        fluxType=problem['flux_type'])
    
    # Initialize the solver
    driver = LinAdvecDriver1D(problem, sp, u0 = 0*sp.x)
    
    # Set the time integrator
    driver.run(FinalTime=problem['tmax'])

    # Logging
    driver.DriverLog(folder_name='cem_3p1')

    # Compute L2 error
    ua = analytical_solution(problem, sp.x)
    L2_error = driver.ComputeL2Error(ua)
    print(f"\n🌐 Erro global na norma L2: {L2_error:.2e}")

    # Plotting the solution
    plot_solution(problem, sp, driver.u_h)


def run_convergence(problem_template):
    """
    Executa um estudo de convergência variando ordem N e número de elementos K.
    Produz tabela estilo LaTeX, com N em linhas, K em colunas e taxa média de convergência.

    Parâmetros
    ----------
    problem_template : dict
        Dicionário base com parâmetros do problema. 'n_order' e 'k_elem' serão sobrescritos.

    Retorno
    -------
    pd.DataFrame
        DataFrame com os dados num formato [N x K] + coluna de taxa.
    """
    N_list = [1, 2, 4, 8]
    K_list = [2, 4, 8, 16, 32, 64]
    table_data = {}

    for N in N_list:
        errors = []
        for K in K_list:
            problem = problem_template.copy()
            # problem['n_order'] = N
            # problem['k_elem'] = K

            sp = DG1D(
                n_order=N,
                mesh=Mesh1D(
                    xmin=0.0,
                    xmax=problem['L'],
                    k_elem=K,
                    boundary_label="Periodic"
                ),
                fluxType=problem['flux_type'])

            driver = LinAdvecDriver1D(problem, sp, u0=0* sp.x)
            driver.run(FinalTime=problem['tmax'])
            ua = analytical_solution(problem, sp.x)
            err = driver.ComputeL2Error(ua)
            errors.append(err)

        # Calcula taxa média de convergência com base nos 3 últimos erros
        rates = []
        for i in range(1, len(errors)):
            if errors[i-1] == 0 or errors[i] == 0:
                rates.append(None)
            else:
                rates.append(np.log(errors[i-1] / errors[i]) / np.log(2))

        # Média das 3 últimas taxas (usando os 3 últimos valores)
        avg_rate = np.mean(rates[-3:]) if all(r is not None for r in rates[-3:]) else None

        # Formata erros com notação científica (1 dígito)
        row = [f"{e:.1E}" for e in errors]
        table_data[N] = row + [f"{avg_rate:.1f}" if avg_rate else "—"]

    # Monta DataFrame
    df = pd.DataFrame.from_dict(table_data, orient='index', columns=[str(k) for k in K_list] + ["Convergence rate"])
    df.index.name = "N\\K"

    # Impressão
    print("\n📊 Tabela Global de Erros L2 e Taxas de Convergência")
    print(df)

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

        u0 = np.sin(problem['kx'] / problem['lmbda'] * sp.x)
        driver = LinAdvecDriver1D(problem, sp, u0)
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
        sp = data['sp']
        x_ana = np.linspace(0, problem['L'], 1000)
        u_ana = analytical_solution(problem, x_ana)

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
        ax.plot(x_ana, u_ana, 'k--', linewidth=1.2, label="Analítica")
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
    PROBLEM = {'name': 'cem_3p18',
        'description': 'Problem 1 - Computational Electromagnetics List 3',
        'flux_type': 'Upwind',  # 'Upwind' or 'Centered'
        'advec_speed': 2*np.pi, # Velocidade de advecção
        'kx': 1,                # Número de onda
        'cfl': 0.75,            # Número de Courant-Friedrichs-Lewy
        'tmax': 10,             # Tempo final da simulação
        'L': 4*np.pi,           # Comprimento do domínio
        'n_order': 2,           # Ordem de interpolação polinomial
        'k_elem': 8             # Número de elementos na malha
    }

    # Test 0 - The solution
    # single_test_solution(PROBLEM)

    # Test 1 - Run convergence study
    run_convergence(PROBLEM)


if __name__ == '__main__':
    main()
