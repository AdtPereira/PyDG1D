#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# pylint: disable=E0611, C0413, C0103, W0401, W0614, W0622
"""
══════════════════════════════════════════════════════════════════════════

ESTE SCRIPT UTILIZA CAMINHOS ABSOLUTOS E NÃO ALTERA O DIRETÓRIO DE TRABALHO

EXECUÇÃO:
cd C:\\git\\PyDG1D
conda activate pyDG1D
python examples\\cem_5\\tmz_cavity.py

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
    'name': 'tmz_cavity',
    'folder': 'cem_5',
    'm': 1,                     # Número de modo
    'n': 1,                     # Número de modo
    'dg': {
        'n_order': 6,               # Ordem de interpolação polinomial
        'flux_type': 'Upwind',      # 'Upwind' or 'Centered'
        'cfl': 1.0,                 # Número de Courant-Friedrichs-Lewy
        'bc': "PEC",                # Condição de contorno: 'PEC', 'SMA'  or 'Periodic
        't0': 0.0,                  # Tempo inicial
        't_final': 10.0,            # Tempo final da simulação
    },
    'domain': {
        'type': 'cubic',            # Tipo de domínio: 'rectangle' ou 'cubic'
        'h': 4.0,                   # Tamanho máximo do elemento da malha
        'Lx': 1.0,                  # Semi-lados do retângulo externo (domínio total)
        'Ly': 1.0,                  # Dimensão total do domínio na direção y
        'Lz': 1.0,                  # Dimensão total do domínio na direção z
        'GID_TFZ': 1,               # Grupo físico para a Total Field Zone (TFZ)
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
    sa = ResonantCavity3D(problem)
    
    # 1. Criar a malha cúbica com dimensão L com Gmsh
    gmsh.initialize()
    mesh_cubeK24(problem)
    EToV = get_EToV(dim=problem['DIM'], index_based=0)
    VX, VY, VZ = extract_VX_VY_VZ(get_nodes_data(dim=problem['DIM']))
    # gmsh.fltk.run()
    gmsh.finalize()

    # 2. Criar o objeto Mesh3D 
    mesh = Mesh3D(vx=VX, vy=VY, vz=VZ, EToV=EToV, boundary_label=problem['dg']['bc'])
    # mesh.plot_mesh(title="mesh_cubeK24.msh", show_vertices=True, alpha=0.15)

    print(f"\nMalha criada com {mesh.number_of_vertices()} vértices e {mesh.number_of_elements()} elementos.")

    # 3. Definir a discretização espacial usando DG3D
    sp = Maxwell3D(n_order=problem['dg']['n_order'], mesh=mesh, fluxType=problem['dg']['flux_type'])
    print(f"\n🔎 Discretização espacial criada com ordem {sp.n_order}, {sp.mesh.number_of_elements()} elementos e {sp.number_of_nodes_per_element()} pontos por elemento.")

    # 4. Criar o driver para resolver a equação de Maxwell
    driver = MaxwellDriver(sp, CFL=problem['dg']['cfl'])

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
            l2_error = sp.compute_L2_error(uh, ua)
            error_data['L2_error'][field_name].append(l2_error)

        driver.step()
        t += driver.dt

    # 4. Plotar resultados
    sa.plot_L2_error(error_data)

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


def main() -> None:
    """Função principal para execução do script."""
    clear_terminal()

    # Test 0 - Validação de um único teste considerando o código MATLAB
    print("\n🔎 Teste de validação com código MATLAB...")
    single_test_validation(PROBLEM)


if __name__ == '__main__':
    main()
    plt.show()
