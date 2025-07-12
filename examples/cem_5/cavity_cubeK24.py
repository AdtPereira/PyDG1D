#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# pylint: disable=E0611, C0413, C0103, W0401, W0614, W0622
"""
══════════════════════════════════════════════════════════════════════════

ESTE SCRIPT UTILIZA CAMINHOS ABSOLUTOS E NÃO ALTERA O DIRETÓRIO DE TRABALHO

EXECUÇÃO:
cd C:\\git\\PyDG1D
conda activate pyDG1D
python examples\\hesthaven\\cavity_cubeK24.py

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

# Adiciona a raiz do projeto ao PYTHONPATH
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..\..')))

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.cm as cm

from maxwell.dg.dg3d import *
from maxwell.dg.dg1d_tools import jacobi_gauss
from maxwell.driver import *
from maxwell.integrators.LSERK4 import *
from maxwell.utils import *
from maxwell.spectralAnalyzer import ResonantCavity3D

from mesher.create_mesh import *
from mesher.plot_mesh import *


MATDATA = loadmat('C:\git\PyDG1D\examplesData\inputs\hesthaven_3D\dados.mat')
BOUNDARY = [{'tag': 101, 'type': 'Dirichlet', 'value': 0.0, 'name': 'Ez_0'}]
MATERIAL = [{'tag': 201, 'name': 'free_space', 'relative_magnetic_permeability': 1, 'relative_electric_permittivity': 1}]   
INFO_GRAPH = {'cell': False, 'nodes': False, 'edges': False, 'edges_numb': False, 'filepath': 'examplesData/inputs/hesthaven/3D_cavity/3D_cavity.svg'}

PROBLEM = {
    'DIM': 3,
    'description': 'Teste do esquema DGTD-Maxwell 3D.',
    'name': 'cavity_cubeK24',
    'folder': 'cem_5',
    'm': 1,                         # Número de modo
    'n': 1,                         # Número de modo
    'dg': {
        'n_order': 3,               # Ordem de interpolação polinomial
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
    
    # 1. Criar a malha retangular com Gmsh
    gmsh.initialize()
    mesh_cubeK24(problem)
    EToV = get_EToV(dim=problem['DIM'], index_based=0)
    VX, VY, VZ = extract_VX_VY_VZ(get_nodes_data(dim=problem['DIM']))
    # gmsh.fltk.run()
    gmsh.finalize()

    # 2. Criar o objeto Mesh3D 
    mesh = Mesh3D(vx=VX, vy=VY, vz=VZ, EToV=EToV, boundary_label=problem['dg']['bc'])
    mesh.plot_mesh(title="mesh_cubeK24.msh", show_vertices=True, alpha=0.15)
    EToE, EToF = mesh.connectivityMatrices()

    print(f"\nMalha criada com {mesh.number_of_vertices()} vértices e {mesh.number_of_elements()} elementos.")
    print(f"\nCoordenadas dos nós (VX, VY, VZ):\n{mesh.vx}\n{mesh.vy}\n{mesh.vz}")
    print(f"\nConectividade dos elementos (EToV):\n{mesh.EToV}")
    print(f"\nConectividade elemento a elemento (EToE):\n{EToE}")
    print(f"\nConectividade elemento a face (EToF):\n{EToF}") 

    # 3. Definir a discretização espacial usando DG3D
    sp = Maxwell3D(n_order=problem['dg']['n_order'], mesh=mesh, fluxType=problem['dg']['flux_type'])
    print(f"\n🔎 Discretização espacial criada com ordem {sp.n_order}, {sp.mesh.number_of_elements()} elementos e {sp.number_of_nodes_per_element()} pontos por elemento.")

    # 4. Coordenadas físicas dos nós (Np x K)
    r, s, t = xyz_to_rst(*set_nodes_in_equilateral_tetrahedron(sp.n_order))
    x, y, z = nodesCoordinates(sp.n_order, mesh)

    # 5. Campo E = (0, 0, Ez)
    Ez = sa.Ez_field(x, y, z, t=sa.t0)
    Ex = np.zeros_like(Ez)
    Ey = np.zeros_like(Ez)

    # 6. Cálculo do rotacional numérico
    xcurlE, ycurlE, zcurlE = Curl3D(
        Ex, Ey, Ez, sp.Dr, sp.Ds, sp.Dt, sp.rx, sp.sx, sp.tx, sp.ry, sp.sy, sp.ty, sp.rz, sp.sz, sp.tz)

    # 7. Cálculo do rotacional analítico
    curlx_analytical = +sa.dyEz_field(x, y, z, t=sa.t0)
    curly_analytical = -sa.dxEz_field(x, y, z, t=sa.t0)

    # 8. Erro quadrático médio
    def rmse(u_num, u_ref):
        assert u_num.shape[0] == sp.number_of_nodes_per_element()
        assert u_num.shape[1] == sp.mesh.number_of_elements()
        return np.sqrt(np.mean((u_num - u_ref)**2))

    rmse_x = rmse(xcurlE, curlx_analytical)
    rmse_y = rmse(ycurlE, curly_analytical)
    L2x_error = sp.compute_L2_error(xcurlE, curlx_analytical)
    L2y_error = sp.compute_L2_error(ycurlE, curly_analytical)
    
    print(f"\n📏 Erro quadrático médio do operador 'Curl3D' (RMSE):")
    print(f"    curl_x  → {rmse_x:.3e}")
    print(f"    curl_y  → {rmse_y:.3e}")
    
    print(f"\n🌐 Erro na norma L2 do operador 'Curl3D':")
    print(f"    curl_x  → {L2x_error:.3e}")
    print(f"    curl_y  → {L2y_error:.3e}")

    # Plot da solução ANALÍTICA
    xi, yi, zi, ezi_dg = sp.interpolate_dg_solution(Ez, resolution=12)
    ezi_analytical = sa.Ez_field(xi, yi, zi, t=np.zeros_like(ezi_dg))
    sa.plot_field(xi, yi, zi, ezi_analytical, title=f"Solução Analítica em t={sa.t_final:.1f}s")


def plot_L2_error_comparison(all_results):
    """
    Gera um gráfico com a evolução temporal do erro L2 para diferentes ordens polinomiais (N),
    comparando os fluxos Upwind e Central.
    """
    fig, axs = plt.subplots(1, 2, figsize=(12, 5), sharey=True)

    # Mapeia os estilos de linha para cada ordem
    styles = {
        4: 'k-.',
        6: 'k:',
        8: 'k-',
        10: 'k--'
    }

    # Loop sobre fluxos
    for ax, flux_type in zip(axs, ['Upwind', 'Centered']):
        flux_data = all_results[flux_type]

        for N, data in sorted(flux_data.items()):
            time = np.array(data['time'])
            l2_errors = np.array(data['L2_error']['Ez'])

            label = f'N={N}'
            linestyle = styles.get(N, 'k--')
            ax.plot(time, l2_errors, linestyle, label=label)

        ax.set_yscale('log')
        ax.set_xlim([0, max(time)])
        ax.set_ylim([1e-12, 1e0])
        ax.set_xlabel('Time')
        ax.set_title(f'{flux_type} flux')

        ax.grid(True, which='both', ls=':')
        ax.legend()

    axs[0].set_ylabel(r'$L^2$ error')

    fig.suptitle(r'Time trace of the discrete $L^2$-error for $E^z$')
    plt.tight_layout(rect=[0, 0, 1, 0.95])


def L2_error_study(problem):
    """
    Executa um estudo de convergência espectral (em N) para dois fluxos (upwind e central)
    aplicados à equação de Maxwell 3D usando DG.

    Retorna:
    --------
    dict contendo os erros L2 e evolução temporal organizados por tipo de fluxo e ordem N.
    """

    # Configurações comuns
    sa = ResonantCavity3D(problem)
    t_final = sa.t_final
    n_list = [4, 6, 8, 10]
    flux_types = ['Upwind', 'Centered']
    all_results = {flux: {} for flux in flux_types}

    # 1. Gerar a malha uma vez
    gmsh.initialize()
    mesh_cubeK6()
    EToV = get_EToV(dim=problem['DIM'], index_based=0)
    VX, VY, VZ = extract_VX_VY_VZ(get_nodes_data(dim=problem['DIM']))
    gmsh.finalize()

    # 2. Loop sobre tipo de fluxo
    for flux_type in flux_types:
        problem['dg']['flux_type'] = flux_type

        # 3. Loop sobre ordens polinomiais
        for n in n_list:
            print(f"\n▶️ Rodando simulação com N={n}, fluxo={flux_type}")

            # Criar malha e solver
            mesh = Mesh3D(vx=VX, vy=VY, vz=VZ, EToV=EToV, boundary_label=problem['dg']['bc'])
            sp = Maxwell3D(n_order=n, mesh=mesh, fluxType=flux_type)
            driver = MaxwellDriver(sp, CFL=problem['dg']['cfl'])

            # Condição inicial
            driver['Ez'][:] = sa.Ez_field(sp.x, sp.y, sp.z, 0)

            # Preparar dados
            analytical_fields = {'Ez': sa.Ez_field}
            error_data = {'time': [], 'L2_error': {key: [] for key in analytical_fields}}
            t = 0.0

            # Loop no tempo
            for _ in range(int(t_final / driver.dt)):
                error_data['time'].append(t)

                for field_name, analytical_fn in analytical_fields.items():
                    uh = driver[field_name]
                    ua = analytical_fn(sp.x, sp.y, sp.z, t)
                    l2_error = sp.compute_L2_error(uh, ua)
                    error_data['L2_error'][field_name].append(l2_error)

                driver.step()
                t += driver.dt

            all_results[flux_type][n] = error_data

    # 4. Plotar resultados
    plot_L2_error_comparison(all_results)
    print("\n🔍 Estudo de convergência espectral concluído.")
    return all_results


def main() -> None:
    """Função principal para execução do script."""
    clear_terminal()

    # Test 0 - Validação de um único teste considerando o código MATLAB
    print("\n🔎 Teste de validação com código MATLAB...")
    single_test_validation(PROBLEM)

    # Test 1 - Convergence study
    print("\n🔎 Estudo de convergência espectral...")
    L2_error_study(PROBLEM)


if __name__ == '__main__':
    main()
    plt.show()
