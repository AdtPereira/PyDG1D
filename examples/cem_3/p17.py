#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# pylint: disable=E0611, C0413, C0103, W0401, W0614, W0622
"""
══════════════════════════════════════════════════════════════════════════

ESTE SCRIPT UTILIZA CAMINHOS ABSOLUTOS E NÃO ALTERA O DIRETÓRIO DE TRABALHO

EXECUÇÃO:
cd C:\\git\\PyDG1D
conda activate pyDG1D
python examples\\cem_3\\p17.py

══════════════════════════════════════════════════════════════════════════

O script está estruturado para garantir que todos os arquivos e pastas usados
sejam acessados por caminhos absolutos, evitando problemas causados por mudanças
do diretório de trabalho (cwd).

PRINCIPAIS CONSIDERAÇÕES SOBRE DIRETÓRIOS:

1️⃣ Diretório raiz do projeto (CWD_ROOT):
    - CWD_ROOT = Path.cwd()
    - Definido automaticamente como o diretório atual quando o script é iniciado.
    - No seu ambiente, o diretório raiz é sempre: C:\\git\\PYDG1D
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
import matplotlib.pyplot as plt

from maxwell.dg.mesh1d import *
from maxwell.dg.dg1d import *
from maxwell.integrators.LSERK4 import *
from maxwell.utils import *
from maxwell.LinAdvecEquation import Driver1D


def gaussian_pulse(problem_data, x, t):
    """
    Define um pulso gaussiano centrado em x0 com largura controlada por beta.
    """
    x0 = problem_data['source']['x0'] * problem_data['L']   # Centro do pulso
    beta = problem_data['source']['beta']                   # Largura do pulso
    a = problem_data['advec_speed']                        # Velocidade de advecção
    L = problem_data['L']                            # Comprimento do domínio
    shifted_x = (x - a * t) % L # Coordenada com deslocamento periódico
    if beta <= 0:
        raise ValueError("O parâmetro beta deve ser positivo para um pulso gaussiano.")
    
    return np.exp(-beta**2 * (shifted_x - x0)**2)


PROBLEM = {
        'description': 'Problem 7 - Computational Electromagnetics List 3',
        'DIM': 1,                  # Dimensão do problema
        'name': 'p17',            # Nome do problema
        'folder': 'cem_3',       # Pasta do problema
        'advec_speed': 0.2,         # Velocidade de advecção
        'flux_type': 'Centered',    # 'Upwind' or 'Centered'
        'cfl': 0.75,                # Número de Courant-Friedrichs-Lewy
        'tmax': 1,                  # Tempo final da simulação
        'kx': 2*np.pi,              # Número de onda
        'lmbda': 1,                 # Comprimento de onda
        'L': 1,                     # Comprimento do domínio
        'n_order': 3,               # Ordem de interpolação polinomial
        'k_elem': 20,               # Número de elementos na malha
        'initial_cond': gaussian_pulse,  # Condição inicial: pulso gaussiano
        'source': {
            'type': 'gaussian_pulse',
            'x0': 0.5,               # Centro do pulso gaussiano
            'beta': 10}              # Largura do pulso gaussiano
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
            k_elem=problem['k_elem']),
        fluxType=problem['flux_type'])
    
    # Initialize the solver
    Driver1D(problem, sp).run_single_test() 


def main() -> None:
    """Função principal para execução do script."""
    clear_terminal()

    # Test 0 - The solution
    single_test_solution(PROBLEM)


if __name__ == '__main__':
    main()
    plt.show()
