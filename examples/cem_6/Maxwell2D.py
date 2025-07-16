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
from maxwell.dg.mesh2d_visualizer import *
from maxwell.dg.mesh2d_creator import *

from maxwell.driver import *
from maxwell.utils import *
from maxwell.integrators.LSERK4 import *

PROBLEM = {
    'DIM': 2,
    'description': 'Teste de convergência do esquema DGTD bidimensional TMz. Hesthaven, p. 205',
    'name': 'Maxwell2D',
    'folder': 'cem_6',
    'dg': {
        'n_order': 4,               # Ordem de interpolação polinomial
        'flux_type': 'Upwind',      # 'Upwind' or 'Centered'
        'cfl': 1.0,                 # Número de Courant-Friedrichs-Lewy
        'bc': "PEC",                # Condição de contorno: 'PEC', 'SMA'  or 'Periodic
        't_final': 10.0,            # Tempo final da simulação
    },
    'domain': {
        'name': 'SquareK16',        # Nome do domínio
        'type': 'square',           # Tipo de domínio: 'square' ou 'cubic'
        'h': 1.0,                   # Tamanho máximo do elemento da malha
        'Lx': 2.0,                  # Semi-lados do domínio externo na direção x
        'Ly': 2.0,                  # Semi-lados do domínio externo na direção y
    },
    'pml': {
        'type': 'uniaxial', 
        'L': 1.0,                   # Largura da camada da PML
        'pml_order': 2,             # Ordem polinomial da PML
        'R': 1E-4                   # Coeficiente de reflexão na interface da PML
    },
}

# Criação da pasta de saída
OUTPUTS = (Path(__file__).parent.parent.parent / 'examplesData' / 'outputs' / PROBLEM["folder"]).resolve()
DATAFILE_PATH = os.path.join(OUTPUTS, f"{PROBLEM['domain']['name']}_pyData")
OUTPUTS.mkdir(parents=True, exist_ok=True)
os.makedirs(DATAFILE_PATH, exist_ok=True)


def resonant_cavity_ez_field(x, y, t):
    ''' Hesthaven's book p. 205 '''
    m = 1
    n = 1 
    n_pi = np.pi * n
    m_pi = np.pi * m
    w = np.pi * np.sqrt(m**2 + n**2)
    return np.sin(m_pi*x) * np.sin(n_pi*y) * np.cos(w*t)


def resonant_cavity_hx_field(x, y, t):
    ''' Hesthaven's book p. 205 '''
    m = 1
    n = 1 
    n_pi = np.pi * n
    m_pi = np.pi * m
    w = np.pi * np.sqrt(m**2 + n**2)
    return - (n_pi / w) * np.sin(m_pi*x) * np.cos(n_pi*y) * np.sin(w*t)

def resonant_cavity_hy_field(x, y, t):
    ''' Hesthaven's book p. 205 '''
    m = 1
    n = 1 
    n_pi = np.pi * n
    m_pi = np.pi * m
    w = np.pi * np.sqrt(m**2 + n**2)
    return (m_pi / w) * np.cos(m_pi*x) * np.sin(n_pi*y) * np.sin(w*t)


def create_mesh_object(problem: dict, filename: str) -> Mesh2D:    
    try:
        with GmshMeshReader(filename, problem_dim=problem['DIM']) as reader:
            reader.read_global_mesh()
            physical_groups = reader.get_physical_group_data()            
            graphicsFilePath = os.path.join(OUTPUTS, f"{PROBLEM['name']}_{PROBLEM['domain']['name']}")                       
            
            mesh2D = Mesh2D(vx=reader.vx, vy=reader.vy, EToV=reader.EToV,
                            physical_groups=physical_groups,
                            boundary_label=problem['dg']['bc'])
            
            pltMesh = Mesh2DVisualizer(mesh2D.vx, mesh2D.vy, mesh2D.EToV, 
                                       physical_groups=mesh2D.physical_groups,
                                       boundary_label=problem['dg']['bc'],
                                       filePath=graphicsFilePath)
            
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


def save_data_to_file(problem: dict, mesh2D: Mesh2D, dg: Maxwell2D) -> None:
    """Salva os dados da malha em um arquivo .npy."""
    try:

        # mesh2D
        np.save(os.path.join(DATAFILE_PATH, "vx.npy"), mesh2D.vx)
        np.save(os.path.join(DATAFILE_PATH, "vy.npy"), mesh2D.vy)
        np.save(os.path.join(DATAFILE_PATH, "EToV.npy"), mesh2D.EToV)
        EToE, EToF = mesh2D.connectivityMatrices()
        np.save(os.path.join(DATAFILE_PATH, "EToE.npy"), EToE)
        np.save(os.path.join(DATAFILE_PATH, "EToF.npy"), EToF)

        # Mapeamento DG-FEM
        np.save(os.path.join(DATAFILE_PATH, "vmapM.npy"), dg.vmapM)
        np.save(os.path.join(DATAFILE_PATH, "vmapP.npy"), dg.vmapP)
        np.save(os.path.join(DATAFILE_PATH, "vmapB.npy"), dg.vmapB)        
        np.save(os.path.join(DATAFILE_PATH, "mapB.npy"), dg.mapB)

        # print(f"vmapM: {dg.vmapM.shape}: \n{dg.vmapM}")
        # print(f"vmapP: {dg.vmapP.shape}: \n{dg.vmapP}")

        # Mapeamento de interface
        np.save(os.path.join(DATAFILE_PATH, "vmapI_tfz_scatter.npy"), dg.vmapI[102][201])
        np.save(os.path.join(DATAFILE_PATH, "vmapI_tfz_total.npy"), dg.vmapI[102][202])
        # print(f'\n vmapI_tfz_scatter = vmapI[102][201]: {dg.vmapI[102][201]}')

        print(f"📁 Dados salvos com sucesso em {DATAFILE_PATH}.")
        
    
    except Exception as e:
        print(f"❌ Erro ao salvar os dados: {e}")
        import traceback
        traceback.print_exc() # Imprime o traceback completo para depuração


def probe_analysis(problem: dict, dg: Maxwell2D) -> MaxwellDriver:
    """
    Inicializa o solver DG, executa a simulação e avalia os campos em um ponto de
    observação (sonda), comparando com a solução analítica de forma otimizada.
    """
    cfl = problem['dg']['cfl']
    final_time = problem['dg']['t_final']
    driver = MaxwellDriver(dg, timeIntegratorType='LSERK4', CFL=cfl)
    driver['Ez'][:] = resonant_cavity_ez_field(dg.x, dg.y, 0)
    print(f"\n🌐 Driver iniciado com discretização temporal {driver.dt:.2e} s.")

    # --- 1. Definição da Sonda (Probe) ---
    probe_point = (0.5, 0.5)
    distances = np.sqrt((dg.x - probe_point[0])**2 + (dg.y - probe_point[1])**2)
    node_idx, elem_idx = np.unravel_index(np.argmin(distances), dg.x.shape)
    probe_coords = (dg.x[node_idx, elem_idx], dg.y[node_idx, elem_idx])
    print(f"\n✔️ Sonda posicionada em (x,y) = ({probe_coords[0]:.4f}, {probe_coords[1]:.4f})")

    # --- 2. Estruturas para Coleta e Plotagem ---
    # Dicionário para armazenar os dados da sonda de forma organizada
    probe_data = {
        'time': [],
        'Ez': {'dg': [], 'an': []},
        'Hx': {'dg': [], 'an': []},
        'Hy': {'dg': [], 'an': []}
    }

    # Dicionário para mapear chaves de campo para suas funções analíticas
    analytical_func = {
        'Ez': resonant_cavity_ez_field,
        'Hx': resonant_cavity_hx_field,
        'Hy': resonant_cavity_hy_field
    }

    # --- 3. Loop Temporal e Coleta de Dados ---
    timeRange = np.arange(0.0, final_time, driver.dt)
    for t in timeRange:
        driver.step()
        current_time = driver.timeIntegrator.time
        probe_data['time'].append(current_time)
        for field, func in analytical_func.items():
            probe_data[field]['dg'].append(driver[field][node_idx, elem_idx])
            probe_data[field]['an'].append(func(probe_coords[0], probe_coords[1], current_time))

    print(f"\n🌐 Tempo final da simulação: {driver.timeIntegrator.time:.1f} s.")

    # --- 4. Plotagem Otimizada dos Dados da Sonda ---
    fig, axes = plt.subplots(3, 1, figsize=(10, 12), sharex=True)
    fig.suptitle(f'2D PEC Cavity Fields at (x,y)=({probe_coords[0]:.2f}, {probe_coords[1]:.2f})', fontsize=12)

    # Lista de configurações para cada subplot, evitando código repetido
    plot_configs = [
        {'key': 'Ez', 'title': '$E_z$', 'ylabel': '$E_z$ (V/m)'},
        {'key': 'Hx', 'title': '$H_x$', 'ylabel': '$H_x$ (A/m)'},
        {'key': 'Hy', 'title': '$H_y$', 'ylabel': '$H_y$ (A/m)'}
    ]

    # Loop para criar os subplots
    for ax, config in zip(axes, plot_configs):
        ax.plot(probe_data['time'], probe_data[config['key']]['an'], 'b', label='Analítica', linewidth=1.0)
        ax.plot(probe_data['time'], probe_data[config['key']]['dg'], 'r--', label='DG-FEM', linewidth=1.0)
        ax.set_title(config['title'])
        ax.set_ylabel(config['ylabel'])
        ax.grid(True)
        ax.legend(loc='upper right')

    axes[-1].set_xlabel('Tempo (s)') # Adiciona rótulo apenas no último subplot
    plt.tight_layout(rect=[0, 0, 1, 0.96]) 

    return driver


def main() -> None:
    """Função principal para execução do script."""
    clear_terminal()
    
    # 1. Criar a malha e o visualizador
    gmshFileName = os.path.join(OUTPUTS, f"{PROBLEM['name']}_{PROBLEM['domain']['name']}.msh")
    SquareK16Domain(PROBLEM, gmshFileName) # Se precisar criar o arquivo .msh    
    
    mesh2D, pltMesh = create_mesh_object(PROBLEM, gmshFileName)
    if mesh2D is None or pltMesh is None:
        print("Erro ao criar a malha ou o visualizador. Verifique os logs.")
        return
    
    # 2. Definir a discretização espacial usando DG2D
    dg = Maxwell2D(n_order=PROBLEM['dg']['n_order'], mesh=mesh2D, fluxType=PROBLEM['dg']['flux_type'])

    # Plotar os pontos de colocação do solver
    # pltMesh.collocation_points(dg, map_type='M', title='Pontos de Colocação (vmapM)')

    # 3. Salvar os dados da malha em um arquivo .npy
    # save_data_to_file(PROBLEM, mesh2D, dg)

    # 4. Inicializar o solver DG
    driver = probe_analysis(PROBLEM, dg)


if __name__ == '__main__':
    main()
    plt.show()