#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# pylint: disable=E0611, C0413, C0103, W0401, W0614, W0622
"""
══════════════════════════════════════════════════════════════════════════

ESTE SCRIPT UTILIZA CAMINHOS ABSOLUTOS E NÃO ALTERA O DIRETÓRIO DE TRABALHO

EXECUÇÃO:
cd C:\\git\\PyDG1D
conda activate pyDG1D
python examples\\cem_5\\cubeK6.py

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
Data: 12/06/2025
"""

import os
import sys
import pandas as pd
import matplotlib.pyplot as plt

# Adiciona a raiz do projeto ao PYTHONPATH
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..\..')))

from maxwell.dg.dg3d import *
from maxwell.driver import *
from maxwell.utils import *
from maxwell.integrators.LSERK4 import *
from maxwell.spectralAnalyzer import ResonantCavity3D

from mesher.create_mesh import *
from mesher.plot_mesh import *


PROBLEM = {
    'DIM': 3,
    'description': 'Teste de convergência do esquema DGTD tridimensional TMz.',
    'name': 'buildMaps_cubeK6',
    'folder': 'cem_5',
    'bc': "PEC",                # Condição de contorno: 'PEC'  or 'Periodic
    'flux_type': 'Upwind',      # 'Upwind' or 'Centered'
    't0': 0.0,                  # Tempo inicial
    'cfl': 1.0,                 # Número de Courant-Friedrichs-Lewy
    'm': 1,                     # Número de modo
    'n': 1,                     # Número de modo
    'L': 1,                     # Dimensão total do domínio
    'n_order': 2,               # Ordem de interpolação polinomial
}


def single_test_validation(PROBLEM) -> None:
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
    sa = ResonantCavity3D(PROBLEM)
    
    # 1. Criar a malha cúbica com dimensão L com Gmsh
    gmsh.initialize()
    mesh_cubeK6()
    EToV = get_EToV(dim=PROBLEM['DIM'], index_based=0)
    VX, VY, VZ = extract_VX_VY_VZ(get_nodes_data(dim=PROBLEM['DIM']))
    # gmsh.fltk.run()
    gmsh.finalize()

    # 2. Criar o objeto Mesh3D 
    mesh = Mesh3D(vx=VX, vy=VY, vz=VZ, EToV=EToV, boundary_label=PROBLEM['bc'])

    # 3. Plotar a malha de diferentes formas
    print("\n🔎 Plotando a malha cúbica K6...")
    # mesh.plot_mesh(title="mesh_cubeK6.msh", show_vertices=True, alpha=0.15)
    # plot_cubeK6_mesh(VX, VY, VZ, EToV)
    # plot_local_cubeK6_mesh(VX, VY, VZ, EToV)

    print(f"\nMalha criada com {mesh.number_of_vertices()} vértices e {mesh.number_of_elements()} elementos.")

    # 3. Definir a discretização espacial usando DG3D
    sp = Maxwell3D(n_order=PROBLEM['n_order'], mesh=mesh, fluxType=PROBLEM['flux_type'])
    print(f"\n🔎 Discretização espacial criada com ordem {sp.n_order}, {sp.mesh.number_of_elements()} elementos e {sp.number_of_nodes_per_element()} pontos por elemento.")

    # 4. Estruturas de Dados DG-FEM
    print(f"\n vmapM (Dim: {sp.vmapM.shape}): \n", sp.vmapM)
    vmapM_3D = sp.vmapM.reshape((sp.n_fp, sp.n_faces, sp.mesh.number_of_elements()), order='F')
    print_3d_matrices(vmapM_3D, title="vmapM")

    print(f"\n vmapP (Dim: {sp.vmapP.shape}): \n", sp.vmapP)
    vmapP_3D = sp.vmapP.reshape((sp.n_fp, sp.n_faces, sp.mesh.number_of_elements()), order='F')
    print_3d_matrices(vmapP_3D, title="vmapP")

    print(f"\n vmapB (Dim: {sp.vmapB.shape}): \n", sp.vmapB)
    print(f"\n mapB (Dim: {sp.mapB.shape}): \n", sp.mapB)

    # 5. Plotar a malha com os ids dos pontos de colocação
    plot_buildMaps_cubeK6(sp)


def main() -> None:
    """Função principal para execução do script."""
    clear_terminal()

    # Test 0 - Validação de um único teste considerando o código MATLAB
    print("\n🔎 Teste de validação com código MATLAB...")
    single_test_validation(PROBLEM)


if __name__ == '__main__':
    main()
    plt.show()
