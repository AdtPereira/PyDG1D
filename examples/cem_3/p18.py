#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# pylint: disable=E0611, C0413, C0103, W0401, W0614, W0622
"""
══════════════════════════════════════════════════════════════════════════

ESTE SCRIPT UTILIZA CAMINHOS ABSOLUTOS E NÃO ALTERA O DIRETÓRIO DE TRABALHO

EXECUÇÃO:
cd C:\\git\\PyDG1D
conda activate pyDG1D
python examples\\cem_3\\p18.py

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
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..\..')))

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from maxwell.dg.mesh1d import *
from maxwell.dg.dg1d import *
from maxwell.integrators.LSERK4 import *
from maxwell.utils import *
from maxwell.advec_driver import LinAdvecDriver1D


def analytical_solution(problem_data, x, t):
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
    return -np.sin(kx * a * (t - x/a))


PROBLEM = {
        'description': 'Problem 8 - Computational Electromagnetics List 3',
        'DIM': 1,               # Dimensão do problema
        'name': 'p18',          # Nome do problema
        'folder': 'cem_3',      # Pasta do problema
        'advec_speed': 2*np.pi, # Velocidade de advecção
        'flux_type': 'Upwind',  # 'Upwind' or 'Centered'
        'cfl': 0.75,            # Número de Courant-Friedrichs-Lewy
        'kx': 1,                # Número de onda
        'tmax': 10,             # Tempo final da simulação
        'L': 4*np.pi,           # Comprimento do domínio
        'n_order': 2,           # Ordem de interpolação polinomial
        'k_elem': 8,            # Número de elementos na malha
        'initial_cond': analytical_solution,   # Condição inicial 
    }


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
    u_ana = analytical_solution(problem, x_ana, problem['tmax'])
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
    plt.savefig(f"examplesData/outputs/{problem['folder']}/{problem['name']}_1D_linear_advection_equation_{sp.fluxType}_flux.png",
                dpi=300, bbox_inches='tight', pad_inches=0.1)


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
    LinAdvecDriver1D(problem, sp).run_single_test() 


def run_convergence(problem):
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
    K_list = [2, 4, 8, 16, 32]
    table_data = {}

    for N in N_list:
        errors = []
        for K in K_list:
            print(f"Running N = {N}, K = {K}...")
            problem = problem.copy()
            sp = DG1D(
                n_order=N,
                mesh=Mesh1D(
                    xmin=0.0,
                    xmax=problem['L'],
                    k_elem=K,
                    boundary_label="Periodic"
                ),
                fluxType=problem['flux_type'])

            driver = LinAdvecDriver1D(problem, sp)
            driver.run(FinalTime=problem['tmax'])
            ua = analytical_solution(problem, sp.x, problem['tmax'])
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
    k_list = [2, 4, 8, 16, 32]
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
    plt.savefig(f"examplesData/outputs/{problem['folder']}/{problem['name']}_temporal_discrete_energy_evolution.png",
                dpi=300, bbox_inches='tight', pad_inches=0.1)

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
        u_ana = analytical_solution(problem, x_ana, problem['tmax'])

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

    # Test 0 - The solution
    single_test_solution(PROBLEM)

    # Test 1 - Run convergence study
    run_convergence(PROBLEM)


if __name__ == '__main__':
    main()
    plt.show()
