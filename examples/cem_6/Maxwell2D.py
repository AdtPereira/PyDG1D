#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# pylint: disable=E0611, C0413, C0103, W0401, W0614, W0622
"""
══════════════════════════════════════════════════════════════════════════

ESTE SCRIPT UTILIZA CAMINHOS ABSOLUTOS E NÃO ALTERA O DIRETÓRIO DE TRABALHO

EXECUÇÃO:
cd C:\\git\\PyDG1D
conda activate pyDG1D
python examples\\cem_6\\Maxwell2D.py

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
from pathlib import Path
import matplotlib.pyplot as plt

# Adiciona a raiz do projeto ao PYTHONPATH
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..\..')))

from maxwell.dg.dg2d import *
from maxwell.dg.mesh2d import *
from maxwell.dg.mesh2dvisualizer import *
from maxwell.driver import *
from maxwell.utils import *
from maxwell.integrators.LSERK4 import *
from mesher.create2dMesh import *

PROBLEM = {
    'DIM': 2,
    'description': 'Estruturas de dados para o problema de Maxwell 2D',
    'name': 'Maxwell2D',
    'folder': 'cem_6',
    'dg': {
        'n_order': 2,               # Ordem de interpolação polinomial
        'flux_type': 'Upwind',      # 'Upwind' or 'Centered'
        'cfl': 1.0,                 # Número de Courant-Friedrichs-Lewy
        'bc': "PEC",                # Condição de contorno: 'PEC', 'SMA'  or 'Periodic
        't_final': 10.0,            # Tempo final da simulação
    },
    'domain': {
        'name': 'SquareK14',        # Nome do domínio
        'type': 'square',           # Tipo de domínio: 'square' ou 'cubic'
        'h': 4.0,                   # Tamanho máximo do elemento da malha
        'Lx': 1.0,                  # Semi-lados do domínio externo na direção x
        'Ly': 1.0,                  # Semi-lados do domínio externo na direção y
    },
}

# Criação da pasta de saída
OUTPUTS = (Path(__file__).parent.parent.parent / 'examplesData' / 'outputs' / PROBLEM["folder"]).resolve()
OUTPUTS.mkdir(parents=True, exist_ok=True)


def create_mesh_object(problem: dict, filename: str) -> Mesh2D:    
    try:
        with GmshMeshReader(filename, problem_dim=problem['DIM']) as reader:
            # 1. Ler a malha global e os dados por grupo físico.
            reader.read_global_mesh()
            physical_groups = reader.get_physical_group_data()
            
            print(f"\n=== {len(physical_groups.keys())} Grupos Físicos Coletados ===")
            print(f"Grupos Físicos (Gmsh): {gmsh.model.getPhysicalGroups()}")
            print(f"Grupos Físicos (Dict): {physical_groups.keys()}")


            for name, data in physical_groups.items():
                print(f"Grupo '{name}': {data['EToV'].shape[0]} elementos.")

            # 2. Inicializar o visualizador com a malha global
            # Cria os objetos de malha usando os dados lidos do reader
            mesh2D = Mesh2D(vx=reader.vx, vy=reader.vy, EToV=reader.EToV,
                            physical_groups=physical_groups, boundary_label=problem['dg']['bc'])
            pltMesh = Mesh2DVisualizer(mesh2D.vx, mesh2D.vy, mesh2D.EToV, 
                                       physical_groups=mesh2D.physical_groups, boundary_label=problem['dg']['bc'])
            
            # Gráfico 1: Malha completa com numeração
            pltMesh.vertices_and_elements(title=f"Malha Completa - {problem['domain']['name']}")
            
            # Gráfico 2: Mapa de grupos físicos com cores
            pltMesh.gmsh_physical_group_map(physical_groups, title=f"Mapa de Grupos Físicos - {problem['domain']['name']}")
            
            print(f"\n=== Malha criada com {mesh2D.number_of_vertices()} vértices e {mesh2D.number_of_elements()} elementos. ===")
            return mesh2D, pltMesh

    except FileNotFoundError:
        print(f"❌ Erro: O arquivo de malha {filename} não foi encontrado.")
        print("Certifique-se de que o arquivo .msh foi criado corretamente.")
        return
    except Exception as e:
        print(f"❌ Ocorreu um erro inesperado: {e}")
        import traceback
        traceback.print_exc() # Imprime o traceback completo para depuração
        return


def create_data_structure(problem: dict, pltMesh: Mesh2DVisualizer, dg: Maxwell2D) -> None:
    """
    Cria a estrutura de dados para o problema de Maxwell 2D.

    Parâmetros
    ----------
    problem : dict
        Dicionário com os parâmetros do problema.

    Retorna
    -------
    None
    """
    print(f"\n🔎 Discretização espacial criada com ordem {dg.n_order}, {dg.mesh.number_of_elements()} elementos e {dg.number_of_nodes_per_element()} pontos por elemento.")
    print(f"\n🔎 Criando estrutura de dados para o problema: {problem['name']}...")

    print(f"\nPhysical Groups (Dict): {dg.mesh.physical_groups.keys()}")
    for key, value in dg.mesh.physical_groups.items():
        print(f"  - Grupo '{value['name']}' (Tag {key}): {value['vx'].shape[0]} vértices, {value['EToV'].shape[0]} elementos de dimensão {value['dim']}.")

    EToE, EToF = dg.mesh.connectivityMatrices()  
    print(f"\nEToV (Dim: {dg.mesh.EToV.shape}): \n", dg.mesh.EToV)
    print(f"\nEToE (Dim: {EToE.shape}): \n", EToE)
    print(f"\nEToF (Dim: {EToF.shape}): \n", EToF)

    # EToG = dg.mesh.connectivityElementsToPhysicalGroups()
    # FToG = dg.mesh.connectivityFacesToPhysicalGroups()
    print(f"\nEToG (Dim: {dg.mesh.EToG.shape}): \n", dg.mesh.EToG)
    print(f"\nFToG (Dim: {dg.mesh.FToG.shape}): \n", dg.mesh.FToG)

    print(f"\n vmapM (Dim: {dg.vmapM.shape}): \n", dg.vmapM)
    # vmapM_2D = dg.vmapM.reshape((dg.n_fp, dg.n_faces, dg.mesh.number_of_elements()), order='F')
    # display_3d_matrices(vmapM_2D, title="vmapM")

    print(f"\n vmapP (Dim: {dg.vmapP.shape}): \n", dg.vmapP)
    # vmapP_2D = dg.vmapP.reshape((dg.n_fp, dg.n_faces, dg.mesh.number_of_elements()), order='F')
    # display_3d_matrices(vmapP_2D, title="vmapP")

    print(f"\n vmapB (Dim: {dg.vmapB.shape}): \n", dg.vmapB)
    print(f"\n mapB (Dim: {dg.mapB.shape}): \n", dg.mapB)  

    # Plotar os pontos de colocação do solver
    pltMesh.collocation_points(dg, map_type='M', title='Pontos de Colocação (vmapM)')


def main() -> None:
    """Função principal para execução do script."""
    clear_terminal()
    
    # 1. Criar a malha e o visualizador
    gmshFileName = os.path.join(OUTPUTS, f"{PROBLEM['name']}_{PROBLEM['domain']['name']}.msh")
    # createSquareK14(gmshFileName) # Se precisar criar o arquivo .msh
    mesh2D, pltMesh = create_mesh_object(PROBLEM, gmshFileName)
    if mesh2D is None or pltMesh is None:
        print("Erro ao criar a malha ou o visualizador. Verifique os logs.")
        return
    
    # 2. Definir a discretização espacial usando DG2D
    dg = Maxwell2D(n_order=PROBLEM['dg']['n_order'], mesh=mesh2D, fluxType=PROBLEM['dg']['flux_type'])

    # 3. Criar a estrutura de dados para o problema
    create_data_structure(PROBLEM, pltMesh, dg)


if __name__ == '__main__':
    main()
    plt.show()