##!/usr/bin/env python3
# -*- coding: utf-8 -*-
# pylint: disable=E0611, C0413, C0103, W0401, W0614, W0622
"""
══════════════════════════════════════════════════════════════════════════

ESTE SCRIPT UTILIZA CAMINHOS ABSOLUTOS E NÃO ALTERA O DIRETÓRIO DE TRABALHO

EXECUÇÃO:
cd C:\\git\\PyDG1D
conda activate pyDG1D
python examples\\cem_4p5_nodal.py

══════════════════════════════════════════════════════════════════════════

Autor: Adilton Pereira
Data: 05/06/2025
"""

import os
import sys

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.tri import Triangulation
from mpl_toolkits.axes_grid1 import make_axes_locatable
from mpl_toolkits.mplot3d import Axes3D

from maxwell.driver import *
from maxwell.dg.mesh1d import *
from maxwell.dg.dg2d import *
from maxwell.dg.dg2d_tools import *
from maxwell.utils import *


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


def classify_nodal_functions(N):
    """
    Classifica os índices das funções nodais (ℓ_i) em vértices, arestas e interior.

    Retorna:
        dict: {'vertices': [...], 'edges': [...], 'interiors': [...]}
    """
    # Total de nós nodais no triângulo
    Np = (N + 1) * (N + 2) // 2

    # Índices dos nós nos vértices
    vertices = [0, 1, 2]  # Sempre os três primeiros

    # Índices de nós nas arestas
    Fmask, fmask1, fmask2, fmask3 = buildFMask(N)
    edge_nodes = np.unique(np.concatenate([fmask1, fmask2, fmask3]))
    edge_nodes = edge_nodes[~np.isin(edge_nodes, vertices)]  # remover vértices

    # Índices restantes são nós interiores
    all_indices = np.arange(Np)
    interior_nodes = np.setdiff1d(all_indices, np.concatenate([vertices, edge_nodes]))

    return {
        'vertices': list(vertices),
        'edges': list(edge_nodes),
        'interiors': list(interior_nodes)
    }


def plot_2d_nodal(ax, r, s, Z, node_coords, k):
    tpc = ax.tripcolor(r, s, Z, shading='gouraud', cmap='coolwarm')
    ax.plot(*node_coords, 'ko', markersize=3)
    vertex_coords, vertex_labels, offsets = triangle_vertex_info()
    draw_vertex_labels_2d(ax, vertex_coords, vertex_labels, offsets)
    ax.set_aspect('equal')
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_title(f'$\\ell_{{{k}}}$', fontsize=10)
    for spine in ax.spines.values():
        spine.set_visible(False)
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    plt.colorbar(tpc, cax=cax)


def plot_3d_nodal(ax, r, s, Z, tri, node_coords, k):
    ax.plot_trisurf(r, s, Z, triangles=tri.triangles, cmap='coolwarm', edgecolor='none')
    draw_triangle_edges_3d(ax, np.array([[[-1, -1], [1, -1]], [[1, -1], [-1, 1]], [[-1, 1], [-1, -1]]]))
    ax.scatter(*node_coords, np.zeros_like(node_coords[0]), color='k', s=10)
    vertex_coords, vertex_labels, offsets = triangle_vertex_info()
    draw_vertex_labels_3d(ax, vertex_coords, vertex_labels, offsets)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_zticks([])
    ax.view_init(elev=30, azim=45)
    ax.set_title(f'$\\ell_{{{k}}}$', fontsize=10)
    for spine in ax.spines.values():
        spine.set_visible(False)


def plot_nodal_basis(R, S, tri, N):
    r_in, s_in = R, S
    r_nodes, s_nodes = xy_to_rs(*reference_triangle_nodes(N))
    node_coords = (r_nodes, s_nodes)

    Np = (N + 1) * (N + 2) // 2
    ncols = int(np.ceil(np.sqrt(Np)))
    nrows = int(np.ceil(Np / ncols))

    V = vandermonde(N, r_nodes, s_nodes)
    V_inv = np.linalg.inv(V)
    V_eval = vandermonde(N, r_in, s_in)

    fig2d, axs2d = plt.subplots(nrows, ncols, figsize=(3 * ncols, 3 * nrows))
    axs2d = axs2d.flatten()
    fig2d.suptitle(f'Funções nodais (tripcolor) para N = {N}', fontsize=12)

    for k in range(Np):
        coeff = np.zeros(Np)
        coeff[k] = 1.0
        Z = V_eval @ (V_inv @ coeff)
        plot_2d_nodal(axs2d[k], r_in, s_in, Z, node_coords, k)
        axs2d[k].collections[0].set_clim(0.0, 1.0)

    for ax in axs2d[Np:]:
        ax.axis('off')
    
    plt.tight_layout()

    fig3d, axs3d = plt.subplots(nrows, ncols, figsize=(3 * ncols, 3 * nrows), subplot_kw={'projection': '3d'})
    axs3d = axs3d.flatten()
    fig3d.suptitle(f'Funções nodais (trisurf) para N = {N}', fontsize=12)

    for k in range(Np):
        coeff = np.zeros(Np)
        coeff[k] = 1.0
        Z = V_eval @ (V_inv @ coeff)
        plot_3d_nodal(axs3d[k], r_in, s_in, Z, tri, node_coords, k)
        axs3d[k].set_zlim(0.0, 1.0)
    
    for ax in axs3d[Np:]:
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


def plot_nodal_basis_on_edges(N):
    r_nodes, s_nodes = xy_to_rs(*reference_triangle_nodes(N))
    V = vandermonde(N, r_nodes, s_nodes)
    Vinv = np.linalg.inv(V)

    t = np.linspace(0, 1, 200)
    edge_labels = ["Aresta 1 (V1→V2)", "Aresta 2 (V2→V3)", "Aresta 3 (V3→V1)"]

    fig, axs = plt.subplots(1, 3, figsize=(15, 4), sharey=True)
    fig.suptitle(f"Funções nodais sobre as arestas (N = {N})", fontsize=12)

    handles_dict = {}

    for edge_id in range(3):
        r_edge, s_edge = edge_rs(t, edge_id)
        V_edge = vandermonde(N, r_edge, s_edge)
        axs[edge_id].set_title(edge_labels[edge_id])
        axs[edge_id].set_xlabel("$t$")
        axs[edge_id].set_xlim(0, 1)
        axs[edge_id].grid(False)

        for k in range(V.shape[1]):
            coeff = np.zeros(V.shape[1])
            coeff[k] = 1.0
            phi_k = V_edge @ (Vinv @ coeff)

            if np.max(np.abs(phi_k)) > 1e-12:  # evita traços nulos
                label = f"$\\ell_{{{k}}}$"
                line, = axs[edge_id].plot(t, phi_k, label=label)
                if label not in handles_dict:
                    handles_dict[label] = line

    axs[0].set_ylabel("$\\ell_i(t)$")

    # Legenda global
    labels, handles = zip(*sorted(handles_dict.items(), key=lambda x: int(x[0].split('{')[1].split('}')[0])))
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
    print("Exemplo de CEM 4.5 - Funções Nodais em Triângulos")
    N_ORDER = 3
    R, S, tri = generate_reference_grid()
    plot_nodal_basis(R, S, tri, N=N_ORDER)
    plot_nodal_basis_on_edges(N=N_ORDER)
    plt.show()