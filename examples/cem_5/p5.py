##!/usr/bin/env python3
# -*- coding: utf-8 -*-
# pylint: disable=E0611, C0413, C0103, W0401, W0614, W0622
"""
══════════════════════════════════════════════════════════════════════════

ESTE SCRIPT UTILIZA CAMINHOS ABSOLUTOS E NÃO ALTERA O DIRETÓRIO DE TRABALHO

EXECUÇÃO:
cd C:\\git\\PyDG1D
conda activate pyDG1D
python examples\\cem_5\\p5_modal.py

══════════════════════════════════════════════════════════════════════════

Autor: Adilton Pereira
Data: 27/06/2025
"""

import os
import sys

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..\..')))

import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial import Delaunay

from maxwell.driver import *
from maxwell.dg.dg3d_tools import *
from maxwell.utils import *


def modal_index(i: int, j: int, k: int, N: int) -> int:
    assert i >= 0 and j >= 0 and k >= 0
    assert i + j + k <= N

    term1 = (11 + 12*N + 3*N**2) * i / 6
    term2 = (2*N + 3) * j / 2
    term3 = k
    term4 = (2 + N) * (i**2) / 2
    term5 = i * j
    term6 = j**2 / 2
    term7 = i**3 / 6

    m = 1 + term1 + term2 + term3 - term4 - term5 - term6 + term7
    return int(np.ceil(m))


def generate_reference_grid_3d(graph_mode=True, N_grid=8):
    # Cria malha cúbica em [-1, 1]^3
    r = np.linspace(-1, 1, N_grid)
    s = np.linspace(-1, 1, N_grid)
    t = np.linspace(-1, 1, N_grid)

    R, S, T = np.meshgrid(r, s, t) 
    R_flat, S_flat, T_flat = R.ravel(), S.ravel(), T.ravel()

    # Máscara: pontos dentro do tetraedro de referência
    mask = (R_flat + S_flat + T_flat <= -1)
    R_in = R_flat[mask]
    S_in = S_flat[mask]
    T_in = T_flat[mask]

    # Coordenadas dos pontos internos
    points = np.vstack((R_in, S_in, T_in)).T

    # Tetraedrização dos pontos internos (para visualização 3D)
    tri = Delaunay(points)

    # Cria triangulação 3D
    if graph_mode:
        fig = plt.figure(figsize=(8, 6))
        ax = fig.add_subplot(111, projection='3d')

        # Pontos internos
        ax.scatter(R_in, S_in, T_in, s=1, alpha=0.6)
        ax.set_title("Pontos no tetraedro de referência")
        ax.set_xlabel("R")
        ax.set_ylabel("S")
        ax.set_zlabel("T")
        ax.set_xlim([-1, 1])
        ax.set_ylim([-1, 1])    
        ax.set_zlim([-1, 1])
        ax.set_box_aspect([1, 1, 1])  # Aspecto igual para os eixos
        ax.set_xticks([-1, 1])
        ax.set_yticks([-1, 1])   
        ax.set_zticks([-1, 1])
        ax.grid(False)

        # Vértices do tetraedro de referência
        v1 = [-1, -1, -1]  # vértice "inferior"
        v2 = [ 1, -1, -1]  # r máximo
        v3 = [-1,  1, -1]  # s máximo
        v4 = [-1, -1,  1]  # t máximo
        vertices = [v1, v2, v3, v4]

        # Arestas do tetraedro: combinações de vértices
        edges = [
            [v1, v2],
            [v1, v3],
            [v1, v4],
            [v2, v3],
            [v2, v4],
            [v3, v4],
        ]

        # Desenha as arestas
        for edge in edges:
            xs, ys, zs = zip(*edge)
            ax.plot(xs, ys, zs, color='k', linestyle='--', linewidth=1.5)

        ax.view_init(elev=30, azim=45)
        plt.tight_layout()

    return R_in, S_in, T_in, tri


def tetrahedron_vertex_info():
    """
    Retorna informações dos vértices do tetraedro de referência:
    - coordenadas (r, s, t)
    - rótulos
    - deslocamentos para plotagem
    """
    vertex_labels = ["V1", "V2", "V3", "V4"]
    vertex_coords = np.array([
        [-1, -1, -1],  # V1
        [ 1, -1, -1],  # V2
        [-1,  1, -1],  # V3
        [-1, -1,  1]   # V4
    ])
    offsets = np.array([
        [-0.15, -0.15, -0.15],  # para V1
        [+0.10, -0.10, -0.10],  # para V2
        [-0.10, +0.10, -0.10],  # para V3
        [-0.10, -0.10, +0.10]   # para V4
    ])
    return vertex_coords, vertex_labels, offsets


def face_points(face_id, Nf):
    t_vals = np.linspace(-1, 1, Nf)
    s_vals = np.linspace(-1, 1, Nf)
    R, S = np.meshgrid(t_vals, s_vals)
    R, S = R.flatten(), S.flatten()

    if face_id == 0:   # Face (r + s + t = -1) → plana superior
        r = R
        s = S
        t = -1 - r - s
    elif face_id == 1:  # Face r = -1
        r = -1 * np.ones_like(R)
        s = R
        t = S
    elif face_id == 2:  # Face s = -1
        r = R
        s = -1 * np.ones_like(R)
        t = S
    elif face_id == 3:  # Face t = -1
        r = R
        s = S
        t = -1 * np.ones_like(R)
    else:
        raise ValueError("Face inválida")

    # Filtra pontos válidos
    mask = (r >= -1) & (s >= -1) & (t >= -1) & (r + s + t <= -1)
    return r[mask], s[mask], t[mask]


def draw_tetrahedron_edges(ax):
    from mpl_toolkits.mplot3d.art3d import Line3DCollection
    vertex_coords, _, _ = tetrahedron_vertex_info()
    edges = [
        [vertex_coords[0], vertex_coords[1]],
        [vertex_coords[0], vertex_coords[2]],
        [vertex_coords[0], vertex_coords[3]],
        [vertex_coords[1], vertex_coords[2]],
        [vertex_coords[1], vertex_coords[3]],
        [vertex_coords[2], vertex_coords[3]],
    ]
    edge_lines = Line3DCollection(edges, colors='k', linewidths=1.0, linestyles='--')
    ax.add_collection3d(edge_lines)


def draw_edge_labels(ax):
    vertex_coords, _, _ = tetrahedron_vertex_info()
    center = np.mean(vertex_coords, axis=0)
    edges = [
        (vertex_coords[0], vertex_coords[1], "e₁"),
        (vertex_coords[0], vertex_coords[2], "e₂"),
        (vertex_coords[0], vertex_coords[3], "e₃"),
        (vertex_coords[1], vertex_coords[2], "e₄"),
        (vertex_coords[1], vertex_coords[3], "e₅"),
        (vertex_coords[2], vertex_coords[3], "e₆"),
    ]
    for v0, v1, label in edges:
        midpoint = 0.5 * (np.array(v0) + np.array(v1))
        e_vec = np.array(v1) - np.array(v0)
        to_outside = midpoint - center
        proj = np.dot(to_outside, e_vec) / np.dot(e_vec, e_vec) * e_vec
        offset = to_outside - proj
        offset = 0.12 * offset / (np.linalg.norm(offset) + 1e-12)
        ax.text(*(midpoint + offset), label, fontsize=8, color='black',
                ha='center', va='center')


def draw_vertex_labels(ax):
    vertex_coords, vertex_labels, vertex_offsets = tetrahedron_vertex_info()
    for (vx, vy, vz), label, offset in zip(vertex_coords, vertex_labels, vertex_offsets):
        ax.text(vx + offset[0], vy + offset[1], vz + offset[2], label,
                fontsize=9, color='black', ha='center', va='center')


def plot_modal_basis_3d(N, N_face=12):
    Np = int((N + 1) * (N + 2) * (N + 3) // 6)
    ncols = int(np.ceil(np.sqrt(Np)))
    nrows = int(np.ceil(Np / ncols))

    fig = plt.figure(figsize=(4 * ncols, 3.5 * nrows))
    fig.suptitle(f'Funções modais nas faces do tetraedro de referência para N = {N}', fontsize=14)

    idx = 0
    for i in range(N + 1):
        for j in range(N + 1 - i):
            for k in range(N + 1 - i - j):
                ax = fig.add_subplot(nrows, ncols, idx + 1, projection='3d')

                draw_tetrahedron_edges(ax)
                draw_edge_labels(ax)
                draw_vertex_labels(ax)

                for face_id in range(4):
                    r, s, t = face_points(face_id, Nf=N_face)
                    a, b, c = rst_to_abc(r, s, t)
                    Z = simplex_polynomial(a, b, c, i, j, k)
                    mask = np.abs(Z) > 1e-3
                    ax.scatter(r[~mask], s[~mask], t[~mask], c='lightgray', s=4, alpha=0.3)
                    ax.scatter(r[mask], s[mask], t[mask], c=Z[mask], cmap='coolwarm', s=8, alpha=0.9)

                index = modal_index(i, j, k, N)
                ax.set_title(f'$\\psi_{{({i},{j},{k})}} = \\psi_{{{index}}}$', fontsize=9)
                ax.set_xlim([-1, 1])
                ax.set_ylim([-1, 1])
                ax.set_zlim([-1, 1])
                ax.set_xticks([])
                ax.set_yticks([])
                ax.set_zticks([])
                ax.view_init(elev=30, azim=35)
                idx += 1

    for k in range(idx, ncols * nrows):
        fig.add_subplot(nrows, ncols, k + 1).axis('off')

    plt.tight_layout()


def classify_nodal_functions_3d(N):
    Np = (N + 1) * (N + 2) * (N + 3) // 6
    vertex_ids = [0, 1, 2, 3]
    edge_nodes = set()
    face_nodes = set()
    all_ids = np.arange(Np)

    sk = 0
    for i in range(N+1):
        for j in range(N+1 - i):
            for k in range(N+1 - i - j):
                l = N - i - j - k
                if i * j * k * l == 0 and (i + j + k + l) == N:
                    if not (i == N or j == N or k == N or l == N):
                        edge_nodes.add(sk)
                    elif sk not in vertex_ids:
                        face_nodes.add(sk)
                sk += 1

    edge_nodes = list(edge_nodes)
    face_nodes = list(face_nodes)
    interior_nodes = list(set(all_ids) - set(vertex_ids) - set(edge_nodes) - set(face_nodes))
    return {
        'vertices': vertex_ids,
        'edges': edge_nodes,
        'faces': face_nodes,
        'interiors': interior_nodes
    }


def create_axes_grid(ncols, nrows, N):
    fig = plt.figure(figsize=(4 * ncols, 3.5 * nrows))
    fig.suptitle(f"Funções nodais no tetraedro para N = {N}", fontsize=14)
    axes = [fig.add_subplot(nrows, ncols, k + 1, projection="3d") for k in range(ncols * nrows)]
    return fig, axes


def evaluate_basis_function(N, k, Vinv, V_eval):
    Np = Vinv.shape[0]
    coeff = np.zeros(Np)
    coeff[k] = 1.0
    return V_eval @ (Vinv @ coeff)


def label_node_type(k, node_types):
    if k in node_types["vertices"]:
        return "(V)"
    elif k in node_types["edges"]:
        return "(E)"
    elif k in node_types["faces"]:
        return "(F)"
    else:
        return "(I)"


def configure_3d_axis(ax, title):
    ax.set_title(title, fontsize=9)
    ax.set_xlim([-1, 1])
    ax.set_ylim([-1, 1])
    ax.set_zlim([-1, 1])
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_zticks([])
    ax.view_init(elev=25, azim=35)


def plot_nodal_basis_3d(N, N_grid=12):
    # Coordenadas nodais e malha de amostragem
    Xn, Yn, Zn = set_nodes_in_equilateral_tetrahedron(N)
    rn, sn, tn = xyz_to_rst(Xn, Yn, Zn)
    R, S, T, _ = generate_reference_grid_3d(graph_mode=False, N_grid=N_grid)

    # Avaliação da base
    V_nodal = vandermonde3D(N, rn, sn, tn)
    Vinv = np.linalg.inv(V_nodal)
    V_eval = vandermonde3D(N, R, S, T)
    node_types = classify_nodal_functions_3d(N)

    # Subplots
    Np = V_nodal.shape[0]
    ncols = int(np.ceil(np.sqrt(Np)))
    nrows = int(np.ceil(Np / ncols))
    fig, axes = create_axes_grid(ncols, nrows, N)

    for k, ax in enumerate(axes):
        if k >= Np:
            ax.axis("off")
            continue

        Zk = evaluate_basis_function(N, k, Vinv, V_eval)
        p = ax.scatter(R, S, T, c=Zk, cmap="coolwarm", s=6, alpha=0.85)

        label = label_node_type(k, node_types)
        configure_3d_axis(ax, fr"$\ell_{{{k}}}$ {label}")

        cb = fig.colorbar(p, ax=ax, shrink=0.6)
        cb.set_ticks([0.0, 0.5, 1.0])
        cb.ax.tick_params(labelsize=6)

    plt.tight_layout()


if __name__ == "__main__":
    clear_terminal()
    print("Exemplo de CEM 4.5 - Funções Modais e Nodais no Tetraedro de Referência")
    R, S, T, tri = generate_reference_grid_3d(graph_mode=False, N_grid=8)    
    plot_modal_basis_3d(N=1)
    plot_nodal_basis_3d(N=1)
    plt.show()