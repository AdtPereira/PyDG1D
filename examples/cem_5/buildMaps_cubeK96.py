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
from maxwell.spectralAnalyzer import UniaxialPml

from mesher.create_mesh import *
from mesher.plot_mesh import *


PROBLEM = {
    'DIM': 3,
    'description': 'Teste da UPML do esquema DGTD tridimensional.',
    'name': 'upml',
    'folder': 'buildMaps_cubeK96',
    'dg': {
        'n_order': 2,               # Ordem de interpolação polinomial
        'flux_type': 'Upwind',      # 'Upwind' or 'Centered'
        'cfl': 1.0,                 # Número de Courant-Friedrichs-Lewy
        'bc': "SMA",                # Condição de contorno: 'PEC', 'SMA'  or 'Periodic
        't_final': 10.0,            # Tempo final da simulação
    },
    'domain': {
        'type': 'cubic',            # Tipo de domínio: 'rectangle' ou 'cubic'
        'h': 4.0,                   # Tamanho máximo do elemento da malha
        'Lx': 2.0,                  # Dimensão total do domínio na direção x
        'Ly': 2.0,                  # Dimensão total do domínio na direção y
        'Lz': 2.0,                  # Dimensão total do domínio na direção z
    },
    'pml': {
        'type': 'uniaxial', 
        'L': 1.0,                   # Largura da camada da PML
        'pml_order': 2,             # Ordem polinomial da PML
        'R': 1E-4                   # Coeficiente de reflexão na interface da PML
    },
    'source': {
        'type': 'gaussian',         # Tipo de fonte: 'gaussian' ou 'derivative'
        'A0': 1.0,                  # Amplitude do pulso gaussiano
        'sigma_r': 0.25,            # Desvio padrão espacial do pulso gaussiano
        'sigma_t': 1.0,             # Desvio padrão temporal do pulso gaussiano
        'x0': (0.0, 0.0, 0.0),      # Posição x do centro do pulso gaussiano
        't0': 0.0                   # Tempo de pico do pulso gaussiano
    }
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
    VX, VY, VZ = extract_VX_VY_VZ(get_nodes_data(dim=problem['DIM']))
    # gmsh.fltk.run()
    gmsh.finalize()

    # 2. Criar o objeto Mesh3D 
    mesh = Mesh3D(vx=VX, vy=VY, vz=VZ, EToV=EToV, boundary_label=problem['dg']['bc'])

    # 3. Criar grupos físicos na malha
    print("\n🔎 Criando grupos físicos na malha cúbica K96...")
    
    # Criar grupo físico Total Field Zone (TFZ)
    GID_PML = 0 # Grupo 'Default'
    GID_TFZ = 1 # Grupo 'TFZ'

    tfzDim = 2 * (problem['domain']['Lx'] - problem['pml']['L'])
    mesh.box_physical_group(
        group_tag=GID_TFZ,
        group_name="TFZ",
        box_dims=(tfzDim, tfzDim, tfzDim),
        center=(0.0, 0.0, 0.0)
    )

    # 3. Relatório final usando os novos métodos de acesso
    print("\n--- Relatório de Grupos Físicos ---")
    print("Grupos definidos:", mesh.group_names)

    for tag, name in mesh.group_names.items():
        # Usa o novo método para buscar os elementos por grupo
        elements_in_group = mesh.get_elements_by_group(tag)
        print(f"  - Grupo '{name}' (Tag {tag}): {len(elements_in_group)} elementos.")

    assigned_tets = len(mesh.get_elements_by_group(GID_TFZ))
    unassigned_tets_count = len(mesh.get_elements_by_group(0))

    assert mesh.number_of_elements() == assigned_tets + unassigned_tets_count, \
        f"Erro: Total de elementos ({mesh.number_of_elements()}) não corresponde à soma dos atribuídos ({assigned_tets}) e não atribuídos ({unassigned_tets_count})."

    # # 3. Plotar a malha de diferentes formas
    # print("\n🔎 Plotando a malha cúbica K96...")
    # mesh.plot_mesh(title="mesh_cubeK96.msh", show_vertices=True, alpha=0.15)

    print(f"\nMalha criada com {mesh.number_of_vertices()} vértices e {mesh.number_of_elements()} elementos.")

    # 3. Definir a discretização espacial usando DG3D
    sp = Maxwell3D(n_order=problem['dg']['n_order'], mesh=mesh, fluxType=problem['dg']['flux_type'])
    print(f"\n🔎 Discretização espacial criada com ordem {sp.n_order}, {sp.mesh.number_of_elements()} elementos e {sp.number_of_nodes_per_element()} pontos por elemento.")

    # 4. Estruturas de Dados DG-FEM
    SHOW_ELEMENTS = [7, 73]  # Elementos a serem exibidos

    print(f"\n🔎 Criando estruturas de dados da malha para os elementos {SHOW_ELEMENTS}...")
    format_matrix(sp.mesh.EToV, title="EToV", elements=SHOW_ELEMENTS)    
    format_matrix(sp.EToE, title="EToE", elements=SHOW_ELEMENTS)
    format_matrix(sp.EToF, title="EToF", elements=SHOW_ELEMENTS)

    print(f"EToG (Dim: {sp.mesh.EToG.shape}): \n", sp.mesh.EToG)

    vmapM_3D = sp.vmapM.reshape((sp.n_fp, sp.n_faces, sp.mesh.number_of_elements()), order='F')
    vmapP_3D = sp.vmapP.reshape((sp.n_fp, sp.n_faces, sp.mesh.number_of_elements()), order='F')
    print_3d_matrices(vmapM_3D, elements=SHOW_ELEMENTS, title="vmapM")
    print_3d_matrices(vmapP_3D, elements=SHOW_ELEMENTS, title="vmapP")

    print(f"\n vmapB (Dim: {sp.vmapB.shape}): \n", sp.vmapB)
    print(f"\n mapB (Dim: {sp.mapB.shape}): \n", sp.mapB)    

    # Supondo uma interface entre os grupos 0 e 1
    sp.buildInterfaceMaps(GID_PML, GID_TFZ)

    # Agora você pode acessar os mapas diretamente
    if hasattr(sp, 'vmapI_G0'):
        print(f"\n🔎 Mapas de interface criados com sucesso entre os grupos {GID_PML} e {GID_TFZ}.")
        
        # Os nós de cada lado
        print(f"\nNós no lado do Grupo 0 (sp.vmapI_G0): {len(sp.vmapI_G0)} nós")
        print(sp.vmapI_G0)
        print(f"\nNós no lado do Grupo 1 (sp.vmapI_G1): {len(sp.vmapI_G1)} nós")
        print(sp.vmapI_G1)

        # Os NOVOS mapas de ÍNDICES
        print(f"\nÍndices em vmapM para o lado G0 (sp.mapI_G0): {len(sp.mapI_G0)} índices")
        print(sp.mapI_G0)
        print(f"\nÍndices em vmapM para o lado G1 (sp.mapI_G1): {len(sp.mapI_G1)} índices")
        print(sp.mapI_G1)

    # Verificação de consistência:
    # Os nós recuperados de vmapM usando mapI_G0 devem ser idênticos a vmapI_G0
    # assert np.array_equal(sp.vmapM[sp.mapI_G0], sp.vmapI_G0)

    # 5. Plotar a malha com os ids dos pontos de colocação
    # Plota a interface da face 0 (+X) do grupo 1 (TFZ)
    plot_buildInterfaceMaps_cubeK96(sp, GROUP_ID=1, No_FACE=0, tfz_dim=tfzDim)
    plot_buildInterfaceMaps_cubeK96(sp, GROUP_ID=0, No_FACE=0, tfz_dim=tfzDim)


def main() -> None:
    """Função principal para execução do script."""
    clear_terminal()

    # Test 0 - Validação de um único teste considerando o código MATLAB
    print("\n🔎 Teste de validação com código MATLAB...")
    single_test_validation(PROBLEM)


if __name__ == '__main__':
    main()
    plt.show()
