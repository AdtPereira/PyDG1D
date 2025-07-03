#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# pylint: disable=E0611, C0413, C0103, W0401, W0614, W0622
"""
══════════════════════════════════════════════════════════════════════════

ESTE SCRIPT UTILIZA CAMINHOS ABSOLUTOS E NÃO ALTERA O DIRETÓRIO DE TRABALHO

EXECUÇÃO:
cd C:\\git\\PyDG1D
conda activate pyDG1D
python examples\\cem_5\\cavity.py

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
from scipy.io import loadmat
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


MATDATA = loadmat('C:\git\PyDG1D\examplesData\inputs\hesthaven_3D\dados.mat')
BOUNDARY = [{'tag': 101, 'type': 'Dirichlet', 'value': 0.0, 'name': 'Ez_0'}]
MATERIAL = [{'tag': 201, 'name': 'free_space', 'relative_magnetic_permeability': 1, 'relative_electric_permittivity': 1}]   
INFO_GRAPH = {'cell': False, 'nodes': False, 'edges': False, 'edges_numb': False, 'filepath': 'examplesData/inputs/hesthaven/3D_cavity/3D_cavity.svg'}
PROBLEM = {
    'DIM': 3,
    'description': 'Teste de convergência do esquema DGTD tridimensional TMz.',
    'name': 'cavity',
    'folder': 'cem_5',
    'bc': "PEC",                # Condição de contorno: 'PEC'  or 'Periodic
    'flux_type': 'Centered',    # 'Upwind' or 'Centered'
    't0': 0.0,                  # Tempo inicial
    'cfl': 1.0,                 # Número de Courant-Friedrichs-Lewy
    'm': 1,                     # Número de modo
    'n': 1,                     # Número de modo
    'L': 1,                     # Dimensão total do domínio
    'n_order': 6,               # Ordem de interpolação polinomial
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
    mesh_cubeK24(L=PROBLEM['L'])
    EToV = get_EToV(dim=PROBLEM['DIM'], index_based=0)
    VX, VY, VZ = extract_VX_VY_VZ(get_nodes_data(dim=PROBLEM['DIM']))
    # gmsh.fltk.run()
    gmsh.finalize()

    # 2. Criar o objeto Mesh3D 
    mesh = Mesh3D(vx=VX, vy=VY, vz=VZ, EToV=EToV, boundary_label=PROBLEM['bc'])
    # mesh.plot_mesh(title="mesh_cubeK24.msh", show_vertices=True, alpha=0.15)

    print(f"\nMalha criada com {mesh.number_of_vertices()} vértices e {mesh.number_of_elements()} elementos.")

    # 3. Definir a discretização espacial usando DG3D
    sp = Maxwell3D(n_order=PROBLEM['n_order'], mesh=mesh, fluxType=PROBLEM['flux_type'])
    print(f"\n🔎 Discretização espacial criada com ordem {sp.n_order}, {sp.mesh.number_of_elements()} elementos e {sp.number_of_nodes_per_element()} pontos por elemento.")

    # 4. Criar o driver para resolver a equação de Maxwell
    driver = MaxwellDriver(sp, CFL=PROBLEM['cfl'])
    
    # Condição inicial
    driver['Ez'][:] = sa.Ez_field(sp.x, sp.y, sp.z, 0)

    # Preparar dados
    analytical_fields = {'Ez': sa.Ez_field}
    error_data = {'time': [], 'L2_error': {key: [] for key in analytical_fields}}
    t = 0.0

    # Loop no tempo
    for _ in range(int(sa.t_final / driver.dt)):
        error_data['time'].append(t)

        for field_name, analytical_fn in analytical_fields.items():
            uh = driver[field_name]
            ua = analytical_fn(sp.x, sp.y, sp.z, t)
            l2_error = compute_L2_error(sp, uh, ua)
            error_data['L2_error'][field_name].append(l2_error)

        driver.step()
        t += driver.dt

    # 4. Plotar resultados
    plot_L2_error(error_data)

    # 5. Interpolação da solução DG
    print("\n🔎 Interpolando a solução final para visualização...")
    xi, yi, zi, ezi = sp.interpolate_dg_solution(driver['Ez'], resolution=12)

    # Calcular a solução analítica nos mesmos pontos e no mesmo tempo
    ezi_analytical = sa.Ez_field(xi, yi, zi, t)
    
    # Plot da solução NUMÉRICA
    sa.plot_field(
        x_i=xi, y_i=yi, z_i=zi, field_data=ezi,
        title=f"Solução Numérica (DG) em t={t:.1f}s"
    )

    # Plot da solução ANALÍTICA
    sa.plot_field(
        x_i=xi, y_i=yi, z_i=zi, field_data=ezi_analytical,
        title=f"Solução Analítica em t={t:.1f}s"
    )

    return driver


def plot_L2_error(error_data) -> None:
    """
    Plota a evolução do erro L2 ao longo do tempo para um ou mais campos.

    A função gera um gráfico com o tempo no eixo x e o erro L2 no eixo y,
    utilizando uma escala logarítmica para melhor visualização de pequenas
    variações de erro.

    Parâmetros
    ----------
    error_data : dict
        Dicionário contendo os dados do erro. Espera-se a seguinte estrutura:
        {
            'time': [t_0, t_1, ...],
            'L2_error': {
                'field_name_1': [error_0, error_1, ...],
                'field_name_2': [error_0, error_1, ...],
                ...
            }
        }
        Onde 'field_name_1' poderia ser 'Ez', por exemplo.

    Retorna
    -------
    None
        A função exibe um gráfico e não retorna nenhum valor.
    """
    time_values = error_data['time']
    l2_error_dict = error_data['L2_error']

    # 1. Criar a figura e os eixos para o plot
    # O uso de subplots é uma boa prática para ter mais controle sobre a figura.
    _, ax = plt.subplots(figsize=(10, 6), dpi=90)

    # 2. Iterar sobre cada campo no dicionário de erros e plotar
    for field_name, errors in l2_error_dict.items():
        # Usamos semilogy para que o eixo y fique em escala logarítmica.
        # Isso é ideal para visualizar ordens de magnitude do erro.
        ax.semilogy(time_values, errors, label=f'Erro L2 - {field_name}', marker='o', markersize=3, linestyle='-')

    # 3. Configurar os detalhes do gráfico para clareza
    ax.set_title('Evolução do Erro $L^2$ ao Longo do Tempo', fontsize=16)
    ax.set_xlabel('Tempo (s)', fontsize=12)
    ax.set_ylabel('Erro $L^2$ (Escala Log)', fontsize=12)
    ax.legend()  # Adiciona a legenda (baseada nos 'labels' definidos no plot)
    ax.grid(True, which='both', linestyle='--', linewidth=0.5) # Adiciona uma grade para facilitar a leitura

    # 4. Ajustar o layout e exibir o gráfico
    plt.tight_layout()
    plt.show()


def main() -> None:
    """Função principal para execução do script."""
    clear_terminal()

    # Test 0 - Validação de um único teste considerando o código MATLAB
    print("\n🔎 Teste de validação com código MATLAB...")
    single_test_validation(PROBLEM)


if __name__ == '__main__':
    main()
    plt.show()
