##!/usr/bin/env python3
# -*- coding: utf-8 -*-
# pylint: disable=E0611, C0413, C0103, W0401, W0614, W0622
"""
══════════════════════════════════════════════════════════════════════════

ESTE SCRIPT UTILIZA CAMINHOS ABSOLUTOS E NÃO ALTERA O DIRETÓRIO DE TRABALHO

EXECUÇÃO:
cd C:\\git\\PyDG1D
conda activate pyDG1D
python examples\\cem_4p5_modal.py

══════════════════════════════════════════════════════════════════════════

Autor: Adilton Pereira
Data: 05/06/2025
"""

import os
import sys

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
    assert i >= 0 and j >= 0 and (i + j) <= N
    return j + (N + 1) * i + 1 - (i * (i - 1)) // 2


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


def triangle_vertex_info():
    vertex_labels = ["V1", "V2", "V3"]
    vertex_coords = np.array([[-1, -1], [1, -1], [-1, 1]])
    offsets = np.array([[-0.05, -0.12], [-0.05, -0.12], [-0.05, 0.12]])
    return vertex_coords, vertex_labels, offsets


def draw_vertex_labels_2d(ax, vertex_coords, vertex_labels, offsets):
    for (vx, vy), offset, label in zip(vertex_coords, offsets, vertex_labels):
        ax.text(vx + offset[0], vy + offset[1], label, fontsize=8, ha='center', va='center')


def draw_vertex_labels_3d(ax, vertex_coords, vertex_labels, offsets):
    for (vx, vy), offset, label in zip(vertex_coords, offsets, vertex_labels):
        ax.text(vx + offset[0], vy + offset[1], 0.0, label, fontsize=8, ha='center', va='center')


def draw_triangle_edges_3d(ax, edges):
    for edge in edges:
        x_edge, y_edge = edge[:, 0], edge[:, 1]
        ax.plot(x_edge, y_edge, np.zeros_like(x_edge), color='k', linewidth=1.5)


def plot_mode_2d(ax, r, s, Z, node_coords, i, j, N):
    tpc = ax.tripcolor(r, s, Z, shading='gouraud', cmap='coolwarm')
    ax.plot(*node_coords, 'ko', markersize=3)
    vertex_coords, vertex_labels, offsets = triangle_vertex_info()
    draw_vertex_labels_2d(ax, vertex_coords, vertex_labels, offsets)
    index = modal_index(i, j, N)
    ax.set_aspect('equal')
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_title(f'$\\psi_{{({i},{j})}} = \\psi_{{{index}}}$', fontsize=10)
    for spine in ax.spines.values():
        spine.set_visible(False)
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    plt.colorbar(tpc, cax=cax)


def plot_mode_3d(ax, r, s, Z, tri, node_coords, i, j, N):
    ax.plot_trisurf(r, s, Z, triangles=tri.triangles, cmap='coolwarm', edgecolor='none')
    draw_triangle_edges_3d(ax, np.array([[[-1, -1], [1, -1]], [[1, -1], [-1, 1]], [[-1, 1], [-1, -1]]]))
    ax.scatter(*node_coords, np.zeros_like(node_coords[0]), color='k', s=10)
    vertex_coords, vertex_labels, offsets = triangle_vertex_info()
    draw_vertex_labels_3d(ax, vertex_coords, vertex_labels, offsets)
    index = modal_index(i, j, N)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_zticks([])
    ax.view_init(elev=30, azim=45)
    ax.set_title(f'$\\psi_{{({i},{j})}} = \\psi_{{{index}}}$', fontsize=10)
    for spine in ax.spines.values():
        spine.set_visible(False)


def plot_modal_basis(R, S, tri, N):
    r_in, s_in = R, S
    a_in, b_in = rs_to_ab(r_in, s_in)
    r_nodes, s_nodes = xy_to_rs(*reference_triangle_nodes(N))
    node_coords = (r_nodes, s_nodes)
    
    Np = (N + 1) * (N + 2) // 2
    ncols = int(np.ceil(np.sqrt(Np)))
    nrows = int(np.ceil(Np / ncols))

    fig2d, axs2d = plt.subplots(nrows, ncols, figsize=(3 * ncols, 3 * nrows))
    fig2d.suptitle(f'Funções modais (tripcolor) for N = {N}', fontsize=12)
    axs2d = axs2d.flatten()

    vmin, vmax = np.inf, -np.inf
    Z_list = []
    for i in range(N + 1):
        for j in range(N - i + 1):
            Z = simplex_polynomial(a_in, b_in, i, j)
            Z_list.append((i, j, Z))
            vmin = min(vmin, Z.min())
            vmax = max(vmax, Z.max())

    for idx, (i, j, Z) in enumerate(Z_list):
        plot_mode_2d(axs2d[idx], r_in, s_in, Z, node_coords, i, j, N)
        axs2d[idx].collections[0].set_clim(vmin, vmax)
        if idx == len(Z_list) - 1:
            divider = make_axes_locatable(axs2d[idx])
            cax = divider.append_axes("right", size="5%", pad=0.05)
            plt.colorbar(axs2d[idx].collections[0], cax=cax)

    for ax in axs2d[len(Z_list):]:
        ax.axis('off')

    plt.tight_layout()

    fig3d, axs3d = plt.subplots(nrows, ncols, figsize=(3 * ncols, 3 * nrows), subplot_kw={'projection': '3d'})
    fig3d.suptitle(f'Funções modais (trisurf) for N = {N}', fontsize=14)

    for idx, (i, j, Z) in enumerate(Z_list):
        ax3d = axs3d.flatten()[idx]
        plot_mode_3d(ax3d, r_in, s_in, Z, tri, node_coords, i, j, N)
        ax3d.set_zlim(vmin, vmax)
    for ax in axs3d.flatten()[len(Z_list):]:
        ax.axis('off')
    plt.tight_layout()


def edge_rs(t, edge_id):
    if edge_id == 0:
        return -1 + 2 * t, -1 * np.ones_like(t)
    elif edge_id == 1:
        return 1 - 2 * t, -1 + 2 * t
    elif edge_id == 2:
        return -1 * np.ones_like(t), 1 - 2 * t
    raise ValueError("ID de aresta inválido")


def modal_basis_on_edges(N, t):
    values_by_edge = []
    global_ymin, global_ymax = np.inf, -np.inf

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

    return values_by_edge, global_ymin, global_ymax


def plot_modal_basis_on_edges(N):
    if N == 0:
        fig, ax = plt.subplots(figsize=(5, 3))
        ax.plot([0, 1], [1 / np.sqrt(2)] * 2, label="$\\psi_{(0,0)} = \\frac{1}{\\sqrt{2}}$")
        ax.set_title("Modo constante no triângulo (N = 0)")
        ax.set_xlabel("Parâmetro $t$ ao longo da aresta")
        ax.set_ylim(0.4, 0.8)
        ax.grid(True)
        ax.legend()
        plt.show()
        return

    t = np.linspace(0, 1, 200)
    edge_labels = ["Aresta 1 (V1→V2)", "Aresta 2 (V2→V3)", "Aresta 3 (V3→V1)"]
    values_by_edge, ymin, ymax = modal_basis_on_edges(N, t)

    fig, axs = plt.subplots(1, 3, figsize=(15, 4), sharey=True)
    fig.suptitle(f"Funções modais ortonormais sobre as arestas (N = {N})", fontsize=12)

    handles_dict = {}
    for edge_id, ax in enumerate(axs):
        for i, j, phi in values_by_edge[edge_id]:
            index = modal_index(i, j, N)
            label = f"$\\psi_{{{index}}}$"
            line, = ax.plot(t, phi, label=label)
            # Armazena apenas se ainda não existe (para evitar duplicatas)
            if label not in handles_dict:
                handles_dict[label] = line       

        ax.set_title(edge_labels[edge_id])
        ax.set_xlabel("$t$")
        ax.set_xlim(0, 1)
        ax.set_ylim(ymin, ymax)
        ax.grid(True)

    axs[0].set_ylabel("$\\psi(t)$")

    # Separar em listas ordenadas
    labels, handles = zip(*sorted(handles_dict.items()))

    # Criar legenda horizontal abaixo do título
    fig.legend(
        handles, labels,
        loc='upper center',
        bbox_to_anchor=(0.5, 0.92),
        ncol=len(labels),
        fontsize=9,
        frameon=True
    )

    plt.tight_layout(rect=[0, 0.05, 1, 0.95])
    plt.show()

if __name__ == "__main__":
    clear_terminal()
    print("Exemplo de CEM 4.5 - Funções Modais em Triângulos")
    N_ORDER = 2
    R, S, tri = generate_reference_grid()
    plot_modal_basis(R, S, tri, N=N_ORDER)
    plot_modal_basis_on_edges(N=N_ORDER)
    plt.show()