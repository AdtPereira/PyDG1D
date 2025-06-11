#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# pylint: disable=E0611, C0413, C0103, W0401, W0614, W0622
"""
══════════════════════════════════════════════════════════════════════════

ESTE SCRIPT UTILIZA CAMINHOS ABSOLUTOS E NÃO ALTERA O DIRETÓRIO DE TRABALHO

EXECUÇÃO:
cd C:\\git\\PyDG1D
conda activate pyDG1D
python examples\\cem_4p6.py

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
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as mtri

from maxwell.dg.dg2d import *
from maxwell.driver import *
from maxwell.integrators.LSERK4 import *
from mesher.create_mesh import *
from maxwell.utils import *

BOUNDARY = [{'tag': 101,
                 'type': 'Dirichlet',
                 'value': 0.0,
                 'name': 'Ez_0'}]
       
MATERIAL = [{'tag': 201,
                 'name': 'free_space',
                 'relative_magnetic_permeability': 1,
                 'relative_electric_permittivity': 1}]   
    

def compute_B_field(N, mesh_data):
    """Calcula o campo B = ∇Az × ẑ a partir de Az = x² + y²"""
    x, y = nodes_coordinates_from_gmsh(N, mesh_data)
    Az = x**2 + y**2
    r, s = xy_to_rs(*set_nodes_in_equilateral_triangle(N))
    Dr, Ds = derivateMatrix(N, r, s)
    rx, sx, ry, sy, _ = geometricFactors(x, y, Dr, Ds)
    dAz_dx, dAz_dy = grad(Dr, Ds, Az, rx, sx, ry, sy)
    Bx = dAz_dy
    By = -dAz_dx
    return x, y, Bx, By


def plot_B_field(x, y, Bx, By, N, mesh_data, scale=30):
    """
    Plota o campo vetorial B sobreposto à malha de Gmsh.

    Parâmetros:
        x, y     : coordenadas dos nós (Np, K)
        Bx, By   : componentes do campo (Np, K)
        N        : ordem do polinômio
        mesh_data: dicionário com 'VX', 'VY', 'EToV'
        scale    : fator de escala para o quiver
    """
    # Malha de fundo
    VX, VY = mesh_data["VX"], mesh_data["VY"]
    EToV = mesh_data["EToV"]

    # Vetores para o gráfico
    X = x.ravel(order='F')
    Y = y.ravel(order='F')
    U = Bx.ravel(order='F')
    V = By.ravel(order='F')

    # Plotar malha com triplot
    _, ax = plt.subplots(figsize=(7, 6))
    triangle = mtri.Triangulation(VX, VY, triangles=EToV)
    ax.triplot(triangle, color="gray", linewidth=0.5)

    # Plotar campo vetorial
    ax.quiver(
        X, Y, U, V,
        angles='xy',
        scale_units='xy',
        scale=scale,
        width=0.003,
        color='blue',
        pivot='middle'
    )

    ax.set_title(fr"$\mathbf{{B}} = \nabla A_z \times \hat{{z}}$ para $N = {N}$")
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.axis("equal")
    ax.grid(False)
    plt.tight_layout()


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
    Bx_exact = 2 * y
    By_exact = -2 * x

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
    r, s = xy_to_rs(*set_nodes_in_equilateral_triangle(N))
    M = mass_matrix(N, r, s)

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

    # PROBLEMA: Definições do problema de advecção linear 1D
    PROBLEM = {'name': 'cem_4p6',
        'folder_name': 'cem_4p6',
        'description': 'Teste de convergência do esquema DGTD bidimensional TMz. Hesthaven, p. 205',
        'json_name': None,
        'flux_type': 'Centered',    # 'Upwind' or 'Centered'
        'cfl': 0.1,                 # Número de Courant-Friedrichs-Lewy
        'n_steps': 40,              # Tempo final da simulação
        'kx': 2*np.pi,              # Número de onda
        'L': 2,                     # Comprimento do domínio
        'n_order': 3                # Ordem de interpolação polinomial
    }
    
    # Criar a malha retangular com Gmsh
    mesh_data = mesh_rectangular_domain(PROBLEM, BOUNDARY, MATERIAL, h=0.5, view_mesh=False, mesh_info=False)

    # Testando com N = 3 e N = 4
    for N in [3, 4]:
        x, y, Bx, By = compute_B_field(N, mesh_data)
        erro_L2 = L2_error_gradient(x, y, Bx, By)
        erro_L2_matrix = L2_error_gradient_mass_matrix(Bx, By, x, y, N)
        print(f"N = {N}: Erro L2 do gradiente: {erro_L2:.2e}. Erro L2 com matriz de massa: {erro_L2_matrix:.2e}")
        plot_B_field(x, y, Bx, By, N, mesh_data)
    plt.show()
    

if __name__ == '__main__':
    main()
