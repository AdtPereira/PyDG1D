#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# pylint: disable=E0611, C0413, C0103, W0401, W0614, W0622
"""
══════════════════════════════════════════════════════════════════════════

ESTE SCRIPT UTILIZA CAMINHOS ABSOLUTOS E NÃO ALTERA O DIRETÓRIO DE TRABALHO

EXECUÇÃO:
cd C:\\git\\PyDG1D
conda activate pyDG1D
python examples\\cem_6\\Maxwell2D_InitialData.py

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
        'n_order': 1,               # Ordem de interpolação polinomial
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


def main() -> None:
    """Função principal para execução do script."""
    clear_terminal()
    gmshFileName = os.path.join(OUTPUTS, f"{PROBLEM['name']}_{PROBLEM['domain']['name']}.msh")
    
    # 1. Criar e/ou ler a malha com Gmsh usando a nova classe GmshMeshReader
    try:
        # Usando o gerenciador de contexto para garantir que o gmsh seja finalizado
        with GmshMeshReader(filename=gmshFileName, problem_dim=PROBLEM['DIM']) as reader:
            # createSquareK14(gmshFileName) # Se precisar criar o arquivo .msh

            # --- Abordagem 1: Ler a malha GLOBAL para o grupo físico 201 ---
            GroupTag = (2, 201)
            print(f"\n--- Lendo a Malha Global para o grupo físico {GroupTag} ---")
            reader.read_global_mesh(GroupTag)

            # --- Abordagem 2: Obter dados por GRUPO FÍSICO ---
            print("\n--- Lendo Dados por Grupo Físico ---")
            print(f"Grupos físicos encontrados: {reader.physical_groups}")
            PhysicalGroups = reader._getPhysicalNodesDict()

            for GroupTag, GroupData in PhysicalGroups.items():
                print(f"\nGrupo Físico {GroupTag}: {GroupData['name']}")
                nodes = GroupData['nodes']
                
                print(f"    Vértices: {len(nodes)}")
                for node, coords in nodes.items():
                    print(f"        - Nó {node}: ({coords[0]}, {coords[1]})")
                
                vx, vy = reader._extractVertices(nodes)
                print(f"\n    Coordenadas extraídas:")
                print(f"    vx={vx}")
                print(f"    vy={vy}")
            
            print("\n")
            print("=" * 60)


            # Exibir os grupos físicos e suas entidades
            print(f"\n--- Entidades Geométricas do Modelo ---")
            print(f"{gmsh.model.getEntities()}")

            print("\n--- Grupos Físicos e suas Entidades ---")
            for GroupDim, GroupTag in reader.physical_groups:
                print(f"\nGrupo Físico {GroupTag}: {gmsh.model.getPhysicalName(GroupDim, GroupTag)}")
                GroupEntitiesTags = gmsh.model.getEntitiesForPhysicalGroup(GroupDim, GroupTag)
                print(f"    Tags das Entidades Geométricas de dimensão {GroupDim}: {GroupEntitiesTags}")

                for EntityTag in GroupEntitiesTags:
                    print(f"    Tag {EntityTag}:")
                    elemTypes, _, elemNodeTags = gmsh.model.mesh.getElements(GroupDim, EntityTag)
                    for elemType, elemNodes in zip(elemTypes, elemNodeTags):
                        elementName, dim, order, numNodes, localNodeCoord, numPrimaryNodes = gmsh.model.mesh.getElementProperties(elemType)
                        k = len(elemNodes) // numNodes                
                        
                        print(f"        Tipo de Elemento da Malha: {elementName}")
                        print(f"        Dimensão: {dim}, Ordem: {order}, Número de Nós: {numNodes}")
                        print(f"        Coordenadas dos Nós Locais: {localNodeCoord}")
                        print(f"        Número de Nós Primários: {numPrimaryNodes}")
                        print(f"        Número de Elementos na malha, k: {k}")
                        print(f"        Nós: {elemNodes}")

                EToV = reader._getEToV((GroupDim, GroupTag))
                print(f"\nEToV (Index0) (Dim: {EToV.shape}): \n{EToV}")

            # 2. Plotar a malha
            # Cria os objetos de malha usando os dados lidos do reader
            mesh2D = Mesh2D(vx=reader.vx, vy=reader.vy, EToV=reader.EToV)
            pltMesh = Mesh2DVisualizer(mesh2D.vx, mesh2D.vy, mesh2D.EToV)
            pltMesh.plot_vertices_and_elements(title=f"{PROBLEM['domain']['name']}")
            print(f"\nMalha criada com {mesh2D.number_of_vertices()} vértices e {mesh2D.number_of_elements()} elementos.")
            print("=" * 60)


    except FileNotFoundError:
        print(f"❌ Erro: O arquivo de malha {gmshFileName} não foi encontrado.")
        print("Certifique-se de que o arquivo .msh foi criado corretamente.")
        return
    except Exception as e:
        print(f"❌ Ocorreu um erro inesperado: {e}")
        return

if __name__ == '__main__':
    main()
    plt.show()
