#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# pylint: disable=E0611, C0413, C0103, W0401, W0614, W0622
"""
══════════════════════════════════════════════════════════════════════════

ESTE SCRIPT UTILIZA CAMINHOS ABSOLUTOS E NÃO ALTERA O DIRETÓRIO DE TRABALHO

EXECUÇÃO:
cd C:\\git\\PyDG1D
conda activate pyDG1D
python examples\\cem_5\\p6.py

══════════════════════════════════════════════════════════════════════════

O script está estruturado para garantir que todos os arquivos e pastas usados
sejam acessados por caminhos absolutos, evitando problemas causados por mudanças
do diretório de trabalho (cwd).

PRINCIPAIS CONSIDERAÇÕES SOBRE DIRETÓRIOS:

1️⃣ Diretório raiz do projeto (CWD_ROOT):
    - CWD_ROOT = Path.cwd()
    - Definido automaticamente como o diretório atual quando o script é iniciado.
    - No seu ambiente, o diretório raiz é sempre:
      C:\\git\\PyDG1D
    - Todas as pastas e arquivos são referenciados em relação a esse diretório.

2️⃣ Estrutura esperada:
    - PyDG1D/
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

from maxwell.dg.dg3d import *
from mesher.create_mesh import *
from mesher.plot_mesh import *
from maxwell.utils import *

BOUNDARY = [{'tag': 101, 'type': 'Dirichlet', 'value': 0.0, 'name': 'Ez_0'}]       
MATERIAL = [{'tag': 201, 'name': 'free_space', 'relative_magnetic_permeability': 1, 'relative_electric_permittivity': 1}]   
PROBLEM = {'DIM': 3, 'name': 'cem_5-p6', 'folder_name': 'cem_5', 'description': '', 'L': 2.0}


def compute_B_field_3d(N, mesh):
    """Calcula o campo B = ∇Az × ẑ a partir de Az = x² + y²"""
    x, y, z = nodes_coordinates(N, mesh)
    Az = x**2 + y**2
    r, s, t = xyz_to_rst(*set_nodes_in_equilateral_tetrahedron(N))
    V = vandermonde3D(N, r, s, t)
    Dr, Ds, Dt = differentiation_matrices_3d(N, r, s, t, V)
    rx, sx, tx, ry, sy, ty, rz, sz, tz, J = geometricFactors3D(x, y, z, Dr, Ds, Dt)
    dAz_dx, dAz_dy, dAz_dz = Grad3D(Az, Dr, Ds, Dt, rx, sx, tx, ry, sy, ty, rz, sz, tz)
    Bx = dAz_dy
    By = -dAz_dx
    return x, y, Bx, By


def L2_error_gradient(x, y, Bx, By):
    """
    Calcula o erro na norma L2 do gradiente numérico B = ∇Az × ẑ,
    comparando com a solução analítica B = (2y, -2x).

    Parâmetros:
        x, y (np.ndarray): coordenadas nodais
        Bx, By (np.ndarray): componentes do campo B numérico

    Retorna:
        erro_L2 (float): norma L2 do erro total
    """
    # Campo exato
    Bx_exact = +2*y
    By_exact = -2*x

    # Erro quadrático em cada componente
    err_Bx = (Bx - Bx_exact)**2
    err_By = (By - By_exact)**2

    # Soma dos erros (norma L2 total)
    err_total = err_Bx + err_By

    # Norma L2: raiz da média dos quadrados dos erros
    erro_L2 = np.sqrt(np.sum(err_total) / err_total.size)

    return erro_L2


def L2_error_gradient_mass_matrix(Bx, By, x, y, N):
    """
    Calcula o erro L2 do gradiente B = ∇Az × ẑ usando matriz de massa.

    Parâmetros:
        Bx, By : np.ndarray
            Componentes numéricas do campo vetorial (shape: [Np, K])
        x, y : np.ndarray
            Coordenadas dos nós (mesmo shape de Bx)
        N : int
            Ordem do polinômio

    Retorna:
        erro_L2 : float
            Norma L2 do erro
    """
    Np, K = x.shape
    r, s, t = xyz_to_rst(*set_nodes_in_equilateral_tetrahedron(N))
    M = mass_matrix(N, r, s, t)

    # Solução exata
    Bx_ex = 2 * y
    By_ex = -2 * x

    # Erro
    eBx = Bx - Bx_ex
    eBy = By - By_ex

    # Integração elemento a elemento
    erro_quad = 0.0
    for k in range(K):
        erro_quad += eBx[:, k] @ M @ eBx[:, k]
        erro_quad += eBy[:, k] @ M @ eBy[:, k]

    return np.sqrt(erro_quad)


def main() -> None:
    """Função principal para execução do script."""
    clear_terminal()

    # 1. Criar a malha retangular com Gmsh
    gmsh.initialize()
    mesh_cubeK6()
    
    # 1.1. Obtém a conectividade dos elementos
    EToV = get_EToV(dim=PROBLEM['DIM'], index_based=0)
    VX, VY, VZ = extract_VX_VY_VZ(get_nodes_data(dim=PROBLEM['DIM']))
    # gmsh.fltk.run()
    gmsh.finalize()

    # 2. Criar o objeto Mesh3D 
    mesh=Mesh3D(vx=VX, vy=VY, vz=VZ, EToV=EToV)
    EToE, EToF = mesh.connectivityMatrices()

    # 3. Plotar a malha
    # plot_tetrahedral_mesh(VX, VY, VZ, EToV)
    # plot_local_tetrahedral_mesh(VX, VY, VZ, EToV)

    print(f"\nMalha criada com {mesh.number_of_vertices()} vértices e {mesh.number_of_elements()} elementos.")
    print(f"\nCoordenadas dos nós (VX, VY, VZ):\n{mesh.vx}\n{mesh.vy}\n{mesh.vz}")
    print(f"\nConectividade dos elementos (EToV):\n{mesh.EToV}")
    print(f"\nConectividade elemento a elemento (EToE):\n{EToE}")
    print(f"\nConectividade elemento a face (EToF):\n{EToF}")

    # Testando com N = 3 e N = 4
    print("\nCalculando o campo B...:")
    for N in [3, 4]:
        x, y, Bx, By = compute_B_field_3d(N, mesh)
        erro_L2 = L2_error_gradient(x, y, Bx, By)
        erro_L2_matrix = L2_error_gradient_mass_matrix(Bx, By, x, y, N)
        print(f"N = {N}: Erro L2 do gradiente: {erro_L2:.2e}. Erro L2 com matriz de massa: {erro_L2_matrix:.2e}")

if __name__ == '__main__':
    main()
    plt.show()
