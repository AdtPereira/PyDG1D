##!/usr/bin/env python3
# -*- coding: utf-8 -*-
# pylint: disable=E0611, C0413, C0103, W0401, W0614, W0622
"""
══════════════════════════════════════════════════════════════════════════

ESTE SCRIPT UTILIZA CAMINHOS ABSOLUTOS E NÃO ALTERA O DIRETÓRIO DE TRABALHO

EXECUÇÃO:
cd C:\\git\\PyDG1D
conda activate pyDG1D
python examples\\cem_4p5.py

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
from matplotlib.tri import Triangulation
from matplotlib import cm
from mpl_toolkits.axes_grid1 import make_axes_locatable
from mpl_toolkits.mplot3d import Axes3D

from maxwell.driver import *
from maxwell.dg.mesh1d import *
from maxwell.dg.dg2d import *
from maxwell.dg.dg2d_tools import *
from maxwell.utils import *


def modal_index(i: int, j: int, N: int) -> int:
    """
    Retorna o índice linear m para a base modal do triângulo de referência,
    dado (i, j) com i, j ≥ 0 e i + j ≤ N.
    
    Parâmetros:
    -----------
    i : int
        Ordem do polinômio em 'a' (primeiro índice).
    j : int
        Ordem do polinômio em 'b' (segundo índice).
    N : int
        Ordem máxima total N = i + j.

    Retorna:
    --------
    m : int
        Índice linear correspondente.
    """
    assert i >= 0 and j >= 0 and (i + j) <= N, "Par (i, j) inválido para a ordem N."
    m = j + (N + 1) * i + 1 - (i * (i - 1)) // 2
    return m


def generate_reference_grid(N_grid=40):
    r = np.linspace(-1, 1, N_grid)
    s = np.linspace(-1, 1, N_grid)
    R, S = np.meshgrid(r, s)
    R_flat, S_flat = R.ravel(), S.ravel()
    mask = (R_flat + S_flat <= 0)
    return R_flat[mask], S_flat[mask], Triangulation(R_flat[mask], S_flat[mask])


def reference_triangle_nodes(N):
    if N == 0:
        return np.array([-1.0]), np.array([-1.0])
    x, y = set_nodes_in_equilateral_triangle(N)
    return x, y


def plot_modal_basis(R, S, tri, N):
    def plot_2d_mode(ax, r_in, s_in, r, s, Z, i, j):
        tpc = ax.tripcolor(r_in, s_in, Z, shading='gouraud', cmap='coolwarm')
        ax.plot(r, s, 'ko', markersize=3)
        
        for (vx, vy), offset, label in zip(vertex_coords, label_offsets_2d, vertex_labels):
            ax.text(vx + offset[0], vy + offset[1], label, fontsize=8, ha='center', va='center')
        
        ax.set_aspect('equal')
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_title(f'$(i,j)=({i},{j})$', fontsize=10)
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="5%", pad=0.05)
        plt.colorbar(tpc, cax=cax)

    def plot_3d_mode(ax, r_in, s_in, Z, r, s, i, j):
        index = modal_index(i, j, N)
        ax.plot_trisurf(r_in, s_in, Z, triangles=tri.triangles, cmap='coolwarm', edgecolor='none')
        
        for edge in edges:
            x_edge, y_edge = edge[:, 0], edge[:, 1]
            ax.plot(x_edge, y_edge, np.zeros_like(x_edge), color='k', linewidth=1.5)
        ax.scatter(r, s, np.zeros_like(r), color='k', s=10)
        
        for (vx, vy), offset, label in zip(vertex_coords, label_offsets_3d, vertex_labels):
            ax.text(vx + offset[0], vy + offset[1], 0.0, label, fontsize=8, ha='center', va='center')        

        ax.text2D(0.05, 0.9, f"$\\psi_{{({i},{j})}} = \\psi_{{{index}}}$", transform=ax.transAxes, fontsize=10)
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_zticks([])
        ax.view_init(elev=30, azim=45)
        ax.set_title(f'$(i,j)=({i},{j})$', fontsize=10)

    edges = np.array([[[-1, -1], [1, -1]], [[1, -1], [-1, 1]], [[-1, 1], [-1, -1]]])
    r_in, s_in = R, S

    vertex_labels = ["(-1, -1)", "(1, -1)", "(-1, 1)"]
    vertex_coords = np.array([[-1, -1], [1, -1], [-1, 1]])
    label_offsets_2d = np.array([[0.0, -0.08], [0.0, -0.08], [-0.05, 0.08]])
    label_offsets_3d = np.array([[0.0, -0.1], [0.0, -0.1], [-0.05, 0.1]])

    x, y = reference_triangle_nodes(N)
    r, s = xy_to_rs(x, y)
    a_in, b_in = rs_to_ab(r_in, s_in)

    # ========= FIGURA 2D ==========
    fig2d = plt.figure(figsize=(3 * (N + 1), 3))
    fig2d.suptitle(f'Funções modais 2D (tripcolor) - N = {N}', fontsize=14)

    for idx, i in enumerate(range(N + 1)):
        j = N - i
        Z = simplex_polynomial(a_in, b_in, i, j)
        ax2d = fig2d.add_subplot(1, N + 1, idx + 1)
        plot_2d_mode(ax2d, r_in, s_in, r, s, Z, i, j)
        for spine in ax2d.spines.values():
            spine.set_visible(False)
    plt.tight_layout()

    # ========= FIGURA 3D ==========
    fig3d = plt.figure(figsize=(3 * (N + 1), 3))
    fig3d.suptitle(f'Funções modais 3D (trisurf) - N = {N}', fontsize=14)

    for idx, i in enumerate(range(N + 1)):
        j = N - i
        Z = simplex_polynomial(a_in, b_in, i, j)
        ax3d = fig3d.add_subplot(1, N + 1, idx + 1, projection='3d')
        plot_3d_mode(ax3d, r_in, s_in, Z, r, s, i, j)
    plt.tight_layout()


def plot_modal_basis_on_edges(N):
    """
    Plota as funções modais ortonormais sobre cada aresta do triângulo
    de referência, usando parametrização geométrica exata das arestas.
    """

    if N == 0:
        fig, ax = plt.subplots(figsize=(5, 3))
        ax.plot([0, 1], [1 / np.sqrt(2)] * 2, label=r"$\psi_{(0,0)} = \frac{1}{\sqrt{2}}$")
        ax.set_title("Modo constante no triângulo (N = 0)")
        ax.set_xlabel("Parâmetro $t$ ao longo da aresta")
        ax.set_ylim(0.4, 0.8)
        ax.grid(True)
        ax.legend()
        plt.tight_layout()
        plt.show()
        return

    # Parametrizações geométricas para cada aresta
    def edge_rs(t, edge_id):
        if edge_id == 0:   # Aresta 1: V1 → V2
            r = -1 + 2 * t
            s = -1 * np.ones_like(t)
        elif edge_id == 1: # Aresta 2: V2 → V3
            r = 1 - 2 * t
            s = -1 + 2 * t
        elif edge_id == 2: # Aresta 3: V3 → V1
            r = -1 * np.ones_like(t)
            s = 1 - 2 * t
        return r, s

    edge_labels = ["Aresta 1 (V1→V2)", "Aresta 2 (V2→V3)", "Aresta 3 (V3→V1)"]

    # Geração de pontos de avaliação (mais densos que os nós)
    t = np.linspace(0, 1, 200)
    global_ymin, global_ymax = np.inf, -np.inf
    values_by_edge = []

    # Avalia todas as funções sobre cada aresta
    for edge_id in range(3):
        r, s = edge_rs(t, edge_id)
        a, b = rs_to_ab(r, s)

        edge_values = []
        for i in range(N + 1):
            for j in range(N - i + 1):
                phi = simplex_polynomial(a, b, i, j)
                edge_values.append((i, j, phi))
                global_ymin = min(global_ymin, phi.min())
                global_ymax = max(global_ymax, phi.max())
        values_by_edge.append(edge_values)

    # Plotagem em subgráficos
    fig, axs = plt.subplots(1, 3, figsize=(15, 4))
    fig.suptitle(f"Funções modais ortonormais sobre as arestas (N = {N})", fontsize=14)

    for edge_id in range(3):
        ax = axs[edge_id]
        for i, j, phi in values_by_edge[edge_id]:
            index = modal_index(i, j, N)
            ax.plot(t, phi, label=rf"$\psi_{{{index}}}$")
        ax.set_title(edge_labels[edge_id])
        ax.set_xlabel("Parâmetro $t$ ao longo da aresta")
        ax.set_ylim(global_ymin, global_ymax)
        ax.grid(True)
        ax.legend(fontsize=7)

    plt.tight_layout()


if __name__ == "__main__":
    clear_terminal()
    print("Exemplo de CEM 4.5 - Funções Modais em Triângulos")

    # Plotar funções modais sobre o triângulo de referência
    N_ORDER = 2
    R, S, tri = generate_reference_grid()
    plot_modal_basis(R, S, tri, N=N_ORDER)

    # Plotar funções modais sobre as arestas do triângulo
    plot_modal_basis_on_edges(N=N_ORDER)
    plt.show()
    
