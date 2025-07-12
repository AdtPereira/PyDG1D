#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# pylint: disable=E0611, C0413, C0103, W0401, W0614, W0622
"""
══════════════════════════════════════════════════════════════════════════

ESTE SCRIPT UTILIZA CAMINHOS ABSOLUTOS E NÃO ALTERA O DIRETÓRIO DE TRABALHO

EXECUÇÃO:
cd C:\\git\\PyDG1D
conda activate pyDG1D
python examples\\cem_5\\cubeK96.py

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

from mesher.create_mesh import *
from mesher.plot_mesh import *


PROBLEM = {
    'DIM': 3,
    'description': 'Teste da UPML do esquema DGTD tridimensional.',
    'name': 'upml',
    'folder': 'cem_5',
    'dg': {
        'n_order': 1,               # Ordem de interpolação polinomial
        'flux_type': 'Upwind',      # 'Upwind' or 'Centered'
        'cfl': 1.0,                 # Número de Courant-Friedrichs-Lewy
        'bc': "SMA",                # Condição de contorno: 'PEC', 'SMA'  or 'Periodic
        't_final': 10.0,            # Tempo final da simulação
    },
    'domain': {
        'type': 'cubic',            # Tipo de domínio: 'rectangle' ou 'cubic'
        'h': 4.0,                   # Tamanho máximo do elemento da malha
        'x0': 1.0,                  # Semi-lados do retângulo intermediário (SFZ)
        'Lx': 2.0,                  # Semi-lados do retângulo externo (domínio total)
        'Ly': 2.0,                  # Dimensão total do domínio na direção y
        'Lz': 3.0,                  # Dimensão total do domínio na direção z
        'GID_TFZ': 1,               # Grupo físico para a Total Field Zone (TFZ)
        'GID_SFZ': 0                # Grupo físico para a Source Field Zone (SFZ)
    },
    'physical_groups': {
        'TFZ': {'name': 'Total Field Zone', 'id': 1},
    },
    'pml': {
        'type': 'uniaxial', 
        'L': 1.0,                   # Largura da camada da PML
        'pml_order': 2,             # Ordem polinomial da PML
        'R': 1E-4                   # Coeficiente de reflexão na interface da PML
    },
}


def single_test_validation(problem) -> None:
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
    # 1. Criar a malha cúbica com dimensão L com Gmsh
    gmsh.initialize()
    mesh_cubeK96_upml(problem)
    EToV = get_EToV(dim=problem['DIM'], index_based=0)
    vx, vy, vz = extract_VX_VY_VZ(get_nodes_data(dim=problem['DIM']))
    #gmsh.fltk.run()
    gmsh.finalize()


    # 2. Criar o objeto Mesh3D 
    mesh = Mesh3D(vx, vy, vz, EToV, boundary_label=problem['dg']['bc'])
    mesh.plot_mesh(title="mesh_cubeK96.msh", show_vertices=True, alpha=0.15)
    print(f"\nMalha criada com {mesh.number_of_vertices()} vértices e {mesh.number_of_elements()} elementos.")


    # 3. Criar grupo físico Total Field Zone (TFZ)
    print("\n🔎 Criando grupos físicos na malha cúbica K96...")
    tfzDim = problem['domain']['x0']        # Semi-lados do retângulo intermediário (TFZ)
    mesh.box_physical_group(
        group_tag=problem['physical_groups']['TFZ']['id'],
        group_name=problem['physical_groups']['TFZ']['name'],
        half_dim=(tfzDim, tfzDim, tfzDim), center=(0.0, 0.0, 0.0), group_info=True)


    # 5. Definir a discretização espacial usando DG3D
    sp = Maxwell3D(n_order=problem['dg']['n_order'], mesh=mesh, fluxType=problem['dg']['flux_type'])
    sp.buildGroupInterfaceMaps()
    print(f"\n🔎 Discretização espacial criada com ordem {sp.n_order}...")
    print(f"{sp.mesh.number_of_elements()} elementos totais, {sp.count_boundary_elements()} elementos de borda e {sp.number_of_nodes_per_element()} pontos por elemento.")
    dg_data_structures(problem, sp)


def dg_data_structures(problem, sp, elements_to_show=[4, 75, 77]) -> None:
    print(f"\n🔎 Criando estruturas de dados da malha para os elementos {elements_to_show}...")    
    display_format_matrix(sp.mesh.EToV, title="EToV", elements=elements_to_show)    
    display_format_matrix(sp.EToE, title="EToE", elements=elements_to_show)
    display_format_matrix(sp.EToF, title="EToF", elements=elements_to_show)
    print(f"\nEToG (Dim: {sp.mesh.EToG.shape}): \n", sp.mesh.EToG)

    vmapM_3D = sp.vmapM.reshape((sp.n_fp, sp.n_faces, sp.mesh.number_of_elements()), order='F')
    vmapP_3D = sp.vmapP.reshape((sp.n_fp, sp.n_faces, sp.mesh.number_of_elements()), order='F')
    vmapI_3D = sp.vmapI.reshape((sp.n_fp, sp.n_faces, sp.mesh.number_of_elements()), order='F')

    display_3d_matrices(vmapM_3D, elements=elements_to_show, title=f"vmapM (Dim: {vmapM_3D.shape})")
    display_3d_matrices(vmapP_3D, elements=elements_to_show, title=f"vmapP (Dim: {vmapP_3D.shape})")
    display_3d_matrices(vmapI_3D, elements=elements_to_show, title=f"vmapI (Dim: {vmapI_3D.shape})")

    print(f"\n vmapB (Dim: {sp.vmapB.shape}): \n", sp.vmapB)
    print(f"\n vmapIM (Dim: {sp.vmapIM.shape}): \n", sp.vmapIM)
    print(f"\n vmapIP (Dim: {sp.vmapIP.shape}): \n", sp.vmapIP)
    print(f"\n mapB (Dim: {sp.mapB.shape}): \n", sp.mapB)   
    print(f"\n mapI (Dim: {sp.mapI.shape}): \n", sp.mapI) 
    
    assert np.array_equal(sp.vmapM[sp.mapI], sp.vmapIM)
    assert np.array_equal(sp.vmapP[sp.mapI], sp.vmapIP)

    print(f"\n vmapIM_G1 (Dim: {sp.vmapIM_G1.shape}): \n", sp.vmapIM_G1)
    print(f"\n vmapIP_G1 (Dim: {sp.vmapIP_G1.shape}): \n", sp.vmapIP_G1)
    print(f"\n vmapIM_G0 (Dim: {sp.vmapIM_G0.shape}): \n", sp.vmapIM_G0)
    print(f"\n vmapIP_G0 (Dim: {sp.vmapIP_G0.shape}): \n", sp.vmapIP_G0)

    plot_interface_elements(sp, group_id=problem['domain']['GID_TFZ'], partner_group_id=0)
    plot_interface_elements(sp, group_id=0, partner_group_id=problem['domain']['GID_TFZ'])


def main() -> None:
    """Função principal para execução do script."""
    clear_terminal()
    single_test_validation(PROBLEM)


if __name__ == '__main__':
    main()
    plt.show()
