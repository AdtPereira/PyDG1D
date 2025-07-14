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
    
    try:
        with GmshMeshReader(filename=gmshFileName, problem_dim=PROBLEM['DIM']) as reader:
            # createSquareK14(gmshFileName) # Descomente se precisar criar o arquivo .msh

            # 1. Ler a malha global (necessário para inicializar o visualizador)
            #    e os dados por grupo físico.
            reader.read_global_mesh()
            group_data = reader.get_physical_group_data()
            
            print("\n--- Dados por Grupo Físico Coletados ---")
            for name, data in group_data.items():
                print(f"Grupo '{name}': {data['EToV'].shape[0]} elementos.")
            print("=" * 60)

            # 2. Inicializar o visualizador com a malha global
            pltMesh = Mesh2DVisualizer(reader.vx, reader.vy, reader.EToV)
            
            # 3. Gerar os dois gráficos
            
            # Gráfico 1: Malha completa com numeração
            pltMesh.plot_vertices_and_elements(title=f"Malha Completa - {PROBLEM['domain']['name']}")
            
            # Gráfico 2: Mapa de grupos físicos com cores
            pltMesh.plot_physical_group_map(group_data, title=f"Mapa de Grupos Físicos - {PROBLEM['domain']['name']}")

    except FileNotFoundError:
        print(f"❌ Erro: O arquivo de malha {gmshFileName} não foi encontrado.")
        print("Certifique-se de que o arquivo .msh foi criado corretamente.")
        return
    except Exception as e:
        print(f"❌ Ocorreu um erro inesperado: {e}")
        import traceback
        traceback.print_exc() # Imprime o traceback completo para depuração
        return

if __name__ == '__main__':
    main()
    plt.show()