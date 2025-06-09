##!/usr/bin/env python3
# -*- coding: utf-8 -*-
# pylint: disable=E0611, C0413, C0103, W0401, W0614, W0622
"""
══════════════════════════════════════════════════════════════════════════

ESTE SCRIPT UTILIZA CAMINHOS ABSOLUTOS E NÃO ALTERA O DIRETÓRIO DE TRABALHO

EXECUÇÃO:
cd C:\\git\\PyDG1D
conda activate pyDG1D
python examples\\cem_4p4.py

══════════════════════════════════════════════════════════════════════════

O script está estruturado para garantir que todos os arquivos e pastas usados
sejam acessados por caminhos absolutos, evitando problemas causados por mudanças
do diretório de trabalho (cwd).

PRINCIPAIS CONSIDERAÇÕES SOBRE DIRETÓRIOS:

1️⃣ Diretório raiz do projeto (CWD_ROOT):
    - CWD_ROOT = Path.cwd()
    - Definido automaticamente como o diretório atual quando o script é iniciado.
    - No seu ambiente, o diretório raiz é sempre:
      C:\\git\\PYDG1D
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
Data: 30/05/2025
"""

import os
import sys

# Adiciona a raiz do projeto ao PYTHONPATH
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

from maxwell.driver import *
from maxwell.dg.mesh1d import *
from maxwell.dg.dg2d import *
from maxwell.dg.dg2d_tools import *
from maxwell.utils import *


def plot_wab_nodes(N_values, json_name_list):
    """ Plota os nós do método Warp & Blend para diferentes valores de N.
    Args:
        N_values (list): Lista de valores de N para os quais os nós serão plotados.
    Returns:
        None: Exibe o gráfico com os nós distribuídos em um triângulo equilátero.
    """
    
    _, axs = plt.subplots(1, len(N_values), figsize=(5 * len(N_values), 5))
    if len(N_values) == 1:
        axs = [axs]

    # Extrai vértices reais a partir de N = 1
    vx, vy = set_nodes_in_equilateral_triangle(1)
    vertices = np.column_stack((vx, vy))
    print(f"Vértices do triângulo equilátero:\n {vertices}")

    for i, N in enumerate(N_values):
        # Extrai vértices de referência
        json_file_name = Path(__file__).parent.parent / 'examplesData' / 'inputs' / 'Nodal_set' / f"{json_name_list[i]}.json"
        vx_ref, vy_ref = extract_webdigitized_data(json_file_name)
        
        x, y = set_nodes_in_equilateral_triangle(N)
        axs[i].scatter(x, y, c='k', s=15, zorder=2)
        axs[i].scatter(
            vx_ref, vy_ref,
            marker='o',
            facecolors='none',      # sem preenchimento
            edgecolors='r',         # cor da borda
            s=40,                   # tamanho do marcador
            zorder=3,
            label='Hesthaven (2008)'
        )

        # Desenha as arestas do triângulo usando os vértices reais
        loop = np.vstack([vertices, vertices[0]])  # fecha o triângulo
        axs[i].plot(loop[:, 0], loop[:, 1], 'b-', lw=1, zorder=1)

        axs[i].set_title(f"Warp & Blend for N = {N}")
        axs[i].set_aspect('equal')
        axs[i].axis('off')

    plt.tight_layout()
    plt.show()


def plot_wab_nodes_reference(N_values):
    """
    Plota os pontos warp and blend mapeados para o triângulo de referência I em coordenadas (r,s).
    """
    _, axs = plt.subplots(1, len(N_values), figsize=(5 * len(N_values), 5))
    if len(N_values) == 1:
        axs = [axs]

    # Triângulo de referência I com vértices (r,s)
    ref_vertices = np.array([[-1, -1], [1, -1], [-1, 1]])
    vertex_labels = ["(-1, -1)", "(1, -1)", "(-1, 1)"]

    for i, N in enumerate(N_values):
        x, y = set_nodes_in_equilateral_triangle(N)
        r, s = xy_to_rs(x, y)  # faz o mapeamento para o triângulo de referência

        axs[i].scatter(r, s, c='k', s=15, zorder=2)

        # desenha borda do triângulo de referência I
        loop = np.vstack([ref_vertices, ref_vertices[0]])  # fecha o triângulo
        axs[i].plot(loop[:, 0], loop[:, 1], 'b-', lw=1.5, zorder=1)

        # anota os vértices com coordenadas
        for (vx, vy), label in zip(ref_vertices, vertex_labels):
            axs[i].annotate(label, xy=(vx, vy), xytext=(vx + 0.05, vy + 0.05), fontsize=10)

        axs[i].set_title(f"Warp & Blend (r,s) - N = {N}")
        axs[i].set_aspect('equal')
        axs[i].axis('off')

    plt.tight_layout()
    plt.show()


def compute_condition_numbers(N_values):
    """ Calcula os números de condição da matriz de Vandermonde para diferentes valores de N.
    Args:
        N_values (list): Lista de valores de N para os quais os números de condição serão calculados.
    Returns:    
        list: Lista de tuplas contendo o valor de N e o número de condição correspondente.
    """
    cond_table = []
    for N in N_values:
        x, y = set_nodes_in_equilateral_triangle(N)
        r, s = xy_to_rs(x, y)
        V = vandermonde(N, r, s)
        cond_number = np.linalg.cond(V)
        cond_table.append((N, cond_number))

    print("\nNúmero de Condição da Matriz de Vandermonde:")
    print("-------------------------------------------")
    print("  N   |  Cond(Vandermonde)")
    print("---------------------------")
    for N, cond in cond_table:
        print(f"  {N:2d}  |  {cond:.2e}") 

    return cond_table


def main() -> None:
    """ Função principal para execução do script. """
    clear_terminal()

    N_plot = [4, 6, 8]

    # Plot distribuição de pontos warp and blend
    plot_wab_nodes(N_plot, json_name_list=['hesthaven_66_a', 'hesthaven_66_b', 'hesthaven_66_c'])

    # Plot pontos warp and blend mapeados para o triângulo de referência I
    plot_wab_nodes_reference(N_plot)

    # Tabela de números de condição da Vandermonde
    compute_condition_numbers([3, 5, 8, 15])    


if __name__ == '__main__':
    main()
