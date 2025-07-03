#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# pylint: disable=E0611, C0413, C0103, W0401, W0614, W0622
"""
══════════════════════════════════════════════════════════════════════════

ESTE SCRIPT UTILIZA CAMINHOS ABSOLUTOS E NÃO ALTERA O DIRETÓRIO DE TRABALHO

EXECUÇÃO:
cd C:\\git\\PyDG1D
conda activate pyDG1D
python examples\\cem_3\\p16.py

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


def resonant_cavity_field(problem, x, t):
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
    kx = problem['kx']
    lamb = problem['lmbda']
    a = problem['advec_speed']
    return np.sin(kx / lamb * (x - a * t))


PROBLEM = {
        'description': 'Problem 6 - Computational Electromagnetics List 3',
        'DIM': 1,                   # Dimensão do problema (1D, 2D, 3D)
        'name': 'p16',              # Nome do problema
        'folder': 'cem_3',          # Pasta de saída para logs e resultados
        'advec_speed': 1,           # Velocidade de advecção
        'flux_type': 'Upwind',      # 'Upwind' or 'Centered'
        'cfl': 0.75,                # Número de Courant-Friedrichs-Lewy
        'tmax': 1E1,                # Tempo final da simulação
        'kx': 2*np.pi,              # Número de onda
        'lmbda': 1,                 # Comprimento de onda
        'L': 1,                     # Comprimento do domínio
        'n_order': 2,               # Ordem de interpolação polinomial
        'k_elem': 8,                # Número de elementos na malha
        'initial_cond': resonant_cavity_field # Função para calcular a condição inicial
    }


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
    LinAdvecDriver1D(problem, sp).run_single_test()    


def run_convergence_simulations(problem: dict, k_list: list):
    """
    Executa simulações para múltiplos valores de K, retornando dados para
    estudo de convergência (erro L2), evolução da energia discreta e soluções finais.

    Parâmetros
    ----------
    problem : dict
        Dicionário com os parâmetros do problema.
    k_list : list of int
        Lista de valores de K (número de elementos).

    Retorno
    -------
    results : list of dict
        Cada dicionário contém:
            'K'         : número de elementos,
            'x'         : coordenadas nodais (Np x K),
            'u_h'       : solução numérica final (Np x K),
            'error_L2'  : erro na norma L2,
            'driver'    : objeto LinAdvecDriver1D com históricos de tempo e energia.
    """
    results = []
    plt.figure(figsize=(10, 6))

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

        # Solução analítica para erro L2
        ua = resonant_cavity_field(problem, sp.x, problem['tmax'])
        error_L2 = driver.ComputeL2Error(ua)

        # Plot energia
        plt.plot(driver.time_history, driver.energy_history, label=fr'$K = {K}$')

        results.append({
            'K': K,
            'x': sp.x.copy(),
            'u_h': driver.uh.copy(),
            'error_L2': error_L2,
            'driver': driver
        })

    plt.title("Evolução temporal da energia discreta para diferentes K")
    plt.xlabel("Tempo")
    plt.ylabel("Energia")
    plt.grid(True)
    plt.legend()
    plt.tight_layout()

    file_name = f"_1D-DG-FEM_EnergyEvolution_{sp.fluxType}_flux.svg"
    full_path = os.path.join(driver.path, f"{problem['name']}{file_name}")
    plt.savefig(full_path, dpi=300, bbox_inches='tight', pad_inches=0.1)
    print(f"📄 Soluções numéricas salvas em: {full_path}")

    # Plotar as soluções numéricas
    file_name = f"_1D-DG-FEM_MultiSolution_{sp.fluxType}_flux.svg"
    driver.plot_multi_solutions(problem, results, file_name)

    return results


def process_convergence_table(results: list):
    errors = [res['error_L2'] for res in results]
    k_list = [res['K'] for res in results]

    rates = [None]
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
    print(f"\n📊 Tabela de Erros e Taxas de Convergência (Ordem N = {PROBLEM['n_order']})")
    print(df.to_string(index=False))
    return df


def main() -> None:
    """Função principal para execução do script."""
    clear_terminal()

    # Teste 0: Solução única
    single_test_solution(PROBLEM)

    # Testes 1 e 2: Convergência, energia e solução final
    k_list = [2, 4, 8, 16, 32]
    results = run_convergence_simulations(PROBLEM, k_list)

    # Convergência
    process_convergence_table(results)   


if __name__ == '__main__':
    main()
    plt.show()
