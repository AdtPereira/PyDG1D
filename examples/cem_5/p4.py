##!/usr/bin/env python3
# -*- coding: utf-8 -*-
# pylint: disable=E0611, C0413, C0103, W0401, W0614, W0622
"""
══════════════════════════════════════════════════════════════════════════

ESTE SCRIPT UTILIZA CAMINHOS ABSOLUTOS E NÃO ALTERA O DIRETÓRIO DE TRABALHO

EXECUÇÃO:
cd C:\\git\\PyDG1D
conda activate pyDG1D
python examples\\cem_5\\p4.py

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
Data: 26/06/2025
"""

import os
import sys

# Adiciona a raiz do projeto ao PYTHONPATH
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..\..')))

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

from maxwell.driver import *
from maxwell.dg.dg3d_tools import *
from maxwell.utils import *


def plot_wab_nodes_3d(N_values):
    """ Plota os nós Warp & Blend no tetraedro equilátero para valores de N. """
    fig = plt.figure(figsize=(5 * len(N_values), 5))

    for i, N in enumerate(N_values):
        X, Y, Z = set_nodes_in_equilateral_tetrahedron(N)

        ax = fig.add_subplot(1, len(N_values), i + 1, projection='3d')
        ax.scatter(X, Y, Z, c='k', s=20, zorder=2)

        # Conecta vértices do tetraedro (opcional: mostrar esqueleto)
        v1 = np.array([-1.0, -1/np.sqrt(3), -1/np.sqrt(6)])
        v2 = np.array([+1.0, -1/np.sqrt(3), -1/np.sqrt(6)])
        v3 = np.array([0.0, +2/np.sqrt(3), -1/np.sqrt(6)])
        v4 = np.array([0.0, 0.0, +3/np.sqrt(6)])
        faces = [[v1, v2, v3], [v1, v2, v4], [v2, v3, v4], [v1, v3, v4]]

        for face in faces:
            loop = np.vstack([face, face[0]])
            ax.plot(loop[:, 0], loop[:, 1], loop[:, 2], 'b-', lw=1)

        ax.set_title(f"Warp & Blend: N = {N}")
        ax.set_xlabel("x")
        ax.set_ylabel("y")
        ax.set_zlabel("z")
        ax.view_init(elev=20, azim=30)
        ax.set_box_aspect([1,1,1])
        ax.axis('off')

    plt.tight_layout()


def identify_faces(r, s, t, tol=1e-4):
    L1 = (1 + t)/2
    L2 = (1 + s)/2
    L3 = -(1 + r + s + t)/2
    L4 = (1 + r)/2
    face_ids = np.full_like(L1, -1, dtype=int)
    face_ids[np.abs(L1) < tol] = 1
    face_ids[np.abs(L2) < tol] = 2
    face_ids[np.abs(L3) < tol] = 3
    face_ids[np.abs(L4) < tol] = 4
    return face_ids


def plot_colored_faces(N):
    X, Y, Z = set_nodes_in_equilateral_tetrahedron(N)
    r, s, t = xyz_to_rst(X, Y, Z)
    face_ids = identify_faces(r, s, t)

    colors = {1: 'red', 2: 'blue', 3: 'green', 4: 'orange'}
    labels = {1: 'Face 1 (L1=0)', 2: 'Face 2 (L2=0)', 3: 'Face 3 (L3=0)', 4: 'Face 4 (L4=0)'}

    fig = plt.figure(figsize=(6, 6))
    ax = fig.add_subplot(111, projection='3d')

    for f in range(1, 5):
        mask = face_ids == f
        ax.scatter(X[mask], Y[mask], Z[mask], label=labels[f], color=colors[f], s=25)

    # nós internos
    interior_mask = face_ids == -1
    ax.scatter(X[interior_mask], Y[interior_mask], Z[interior_mask], color='gray', s=20, label="Interior")

    # vértices do tetraedro
    v1 = np.array([-1.0, -1/np.sqrt(3), -1/np.sqrt(6)])
    v2 = np.array([+1.0, -1/np.sqrt(3), -1/np.sqrt(6)])
    v3 = np.array([0.0, +2/np.sqrt(3), -1/np.sqrt(6)])
    v4 = np.array([0.0, 0.0, +3/np.sqrt(6)])
    vertices = [v1, v2, v3, v4]
    vlabels = ["V1", "V2", "V3", "V4"]

    for i, v in enumerate(vertices):
        ax.text(v[0]+0.1, v[1]+0.1, v[2]+0.1, vlabels[i], fontsize=10, color='black')


    # adiciona arestas do tetraedro
    v1 = np.array([-1.0, -1/np.sqrt(3), -1/np.sqrt(6)])
    v2 = np.array([+1.0, -1/np.sqrt(3), -1/np.sqrt(6)])
    v3 = np.array([ 0.0, +2/np.sqrt(3), -1/np.sqrt(6)])
    v4 = np.array([ 0.0, 0.0, +3/np.sqrt(6)])
    edges = [(v1, v2), (v2, v3), (v3, v1), (v1, v4), (v2, v4), (v3, v4)]

    for edge in edges:
        pts = np.vstack(edge)
        ax.plot(pts[:,0], pts[:,1], pts[:,2], 'k--', lw=1)

    ax.set_title(f"Node coordinates in equilateral tetrahedron (Warp & Blend) por N = {N}")
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_zlabel("z")
    ax.legend()
    ax.grid(False)
    ax.view_init(elev=20, azim=30)
    ax.set_box_aspect([1,1,1])

    # Caixa de texto com coordenadas fora do gráfico
    coord_text = '\n'.join([
        f"V1 = ({v1[0]:.1f}, {v1[1]:.1f}, {v1[2]:.1f})",
        f"V2 = ({v2[0]:.1f}, {v2[1]:.1f}, {v2[2]:.1f})",
        f"V3 = ({v3[0]:.1f}, {v3[1]:.1f}, {v3[2]:.1f})",
        f"V4 = ({v4[0]:.1f}, {v4[1]:.1f}, {v4[2]:.1f})",
    ])
    plt.gcf().text(0.75, 0.25, coord_text, fontsize=10,
                   bbox=dict(facecolor='white', edgecolor='black'))
    
    
    plt.tight_layout()


def main():
    """ Função principal para executar o script. """
    clear_terminal()

    # print("Plotando nós Warp & Blend no tetraedro equilátero...")  
    N_values = [2, 3, 4]
    plot_wab_nodes_3d(N_values)

    print("Plotando faces coloridas com interior mínimo para N = 4...")
    plot_colored_faces(N=4)

    # Ordens a testar
    orders = [3, 5, 8, 15]
    condition_numbers = []

    for N in orders:
        r, s, t = set_nodes_in_equilateral_tetrahedron(N)
        V = vandermonde(N, r, s, t)
        cond = np.linalg.cond(V)
        condition_numbers.append((N, cond))
    
    df = pd.DataFrame(condition_numbers, columns=["Ordem N", "Número de Condição"])
    print("\nNúmero de condição de V3D para diferentes ordens N:")
    print(df)


if __name__ == "__main__":
    main()
    plt.show()
