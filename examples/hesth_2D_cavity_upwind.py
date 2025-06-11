#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# pylint: disable=E0611, C0413, C0103, W0401, W0614, W0622
"""
══════════════════════════════════════════════════════════════════════════

ESTE SCRIPT UTILIZA CAMINHOS ABSOLUTOS E NÃO ALTERA O DIRETÓRIO DE TRABALHO

EXECUÇÃO:
cd C:\\git\\PyDG1D
conda activate pyDG1D
python examples\\hesth_2D_cavity.py

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
Data: 26/05/2025
"""

import os
import sys

# Adiciona a raiz do projeto ao PYTHONPATH
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

import numpy as np
import pandas as pd
from pathlib import Path
import matplotlib.pyplot as plt

from maxwell.dg.dg2d import *
from maxwell.driver import *
from maxwell.integrators.LSERK4 import *
from mesher.create_mesh import *
from maxwell.utils import *


BOUNDARY = [{'tag': 101,
                 'type': 'Dirichlet',
                 'value': 0.0,
                 'name': 'Ez_0'}]
   
    
MATERIAL = [{'tag': 201,
                 'name': 'free_space',
                 'relative_magnetic_permeability': 1,
                 'relative_electric_permittivity': 1}]   
    

def resonant_cavity_ez_field(x, y, t):
    ''' Hesthaven's book p. 205 '''
    m = 1
    n = 1 
    w = np.pi * np.sqrt(m**2 + n**2)
    return np.sin(m*np.pi*x)*np.sin(n*np.pi*y)*np.cos(w*t)


def single_test_solution(problem_data) -> None:
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
    # Criar a malha retangular com Gmsh
    mesh_data = mesh_rectangular_domain(problem_data, BOUNDARY, MATERIAL, h=1, view_mesh=False, mesh_info=False) 

    # Visualizar a malha criada
    INFO_GRAPH = {'cell': False, 'nodes': False, 'edges': False, 'edges_numb': False,
                    'filepath': 'examplesData/inputs/hesth_2D_cavity/hesth_2D_cavity.svg'}
    plot_triangular_mesh(INFO_GRAPH, mesh_data)

    # Definir a discretização espacial usando DG1D
    sp = Maxwell2D(
        n_order=problem_data['n_order'],
        mesh=Mesh2D(
            vx=mesh_data['VX'],
            vy=mesh_data['VY'],
            EToV=mesh_data['EToV'],
            boundary_label="PEC"),
        fluxType=problem_data['flux_type'])
    
    # Initialize the solver
    driver = MaxwellDriver(sp, CFL=problem_data['cfl'])
    driver['Ez'][:] = resonant_cavity_ez_field(sp.x, sp.y, 0)
    
    # Set the time integrator
    driver.run(problem_data['n_steps'] * driver.dt)
    print(f"\n🌐 Tempo final da simulação: {driver.timeIntegrator.time:.2f} s")
    
    ez_expected = resonant_cavity_ez_field(sp.x, sp.y, driver.timeIntegrator.time)
    R = np.corrcoef(ez_expected, driver['Ez'])
    print(f"🌐 Correlação entre solução esperada e solução numérica: {R[0,1]:.2f}")


def compute_L2_error(sp, uh, ua):
    """
    Calcula o erro global na norma L2 usando errᵗ diag(J) M err para cada elemento.

    Parâmetros
    ----------
    sp : DG1D
        Objeto com a discretização espacial.
    u_h : ndarray
        Solução numérica final do método DG (Np x K).
    ua : ndarray
        Solução analítica (Np x K).

    Retorno
    -------
    float
        Erro global na norma L2.
    """
    err = ua - uh
    M = sp.mass 
    K = sp.mesh.number_of_elements()
    errL2_local = np.zeros(K)

    for k in range(K):
        Jk = np.diag(sp.jacobian[:, k])           # (Np x Np)
        ek = err[:, k][:, np.newaxis]             # (Np x 1)
        errL2_local[k] = (ek.T @ Jk @ M @ ek)[0, 0]
    return np.sqrt(np.sum(errL2_local))


def plot_L2_error_over_time(error_data: dict) -> None:
    """
    Plota a evolução do erro L2 ao longo do tempo para cada campo registrado.

    Parâmetros
    ----------
    error_data : dict
        Dicionário contendo 'time' e 'L2_error' com os nomes dos campos.
    """
    plt.figure(figsize=(8, 5))
    
    time = error_data['time']
    for field_name, errors in error_data['L2_error'].items():
        plt.plot(time, errors, label=f'Erro L2 - {field_name}', marker='o')
    
    plt.xlabel('Tempo')
    plt.ylabel('Erro $L^2$')
    plt.title('Evolução do erro $L^2$ ao longo do tempo')
    plt.grid(True, which="both", ls='--', alpha=0.6)
    plt.legend()
    plt.tight_layout()
    plt.show()


def run_L2_error(problem_data):
    """ 
    Executa um estudo de convergência para a equação de advecção linear 1D.
    Parâmetros
    ----------
    N : int
        Ordem de interpolação polinomial para o método DG1D.
    
    Retorno
    -------
    pd.DataFrame
    DataFrame contendo os erros L2 e taxas de convergência para diferentes
    números de elementos (K).
    """
    k_list = [0.1]
    for k in k_list:
        # Atualiza o número de elementos na malha
        mesh_data = mesh_rectangular_domain(problem_data, BOUNDARY, MATERIAL, h=k, auto_save=False)

        # Generate Spatial Discretization
        sp = Maxwell2D(
            n_order=problem_data['n_order'],
            mesh=Mesh2D(
                vx=mesh_data['VX'],
                vy=mesh_data['VY'],
                EToV=mesh_data['EToV'],
                boundary_label="PEC"),
            fluxType=problem_data['flux_type'])
        
        # Initialize the solver
        driver = MaxwellDriver(sp, CFL=problem_data['cfl'])
        driver['Ez'][:] = resonant_cavity_ez_field(sp.x, sp.y, 0)
        
        analytical_fields={'Ez': resonant_cavity_ez_field}
        error_data = {'time': [], 'L2_error': {key: [] for key in analytical_fields}}
        t = 0.0

        for _ in range(int(problem_data['n_steps'])):
            error_data['time'].append(t)
            for field_name, analytical_fn in analytical_fields.items():
                uh = driver[field_name]
                ua = analytical_fn(sp.x, sp.y, t)
                l2_error = compute_L2_error(sp, uh, ua)
                error_data['L2_error'][field_name].append(l2_error)

            driver.step()
            t += driver.dt
    
    plot_L2_error_over_time(error_data)
    return error_data


def compute_error_data(problem_data, N, h):
    """Executa a simulação para (N, K) e retorna o erro L2 máximo."""
    # Atualiza o número de elementos na malha
    mesh_data = mesh_rectangular_domain(problem_data, BOUNDARY, MATERIAL, h, auto_save=False)

    # Generate Spatial Discretization
    sp = Maxwell2D(
        n_order=N,
        mesh=Mesh2D(
            vx=mesh_data['VX'],
            vy=mesh_data['VY'],
            EToV=mesh_data['EToV'],
            boundary_label="PEC"),
        fluxType=problem_data['flux_type'])

    # Initialize the solver
    driver = MaxwellDriver(sp, CFL=problem_data['cfl'])
    driver['Ez'][:] = resonant_cavity_ez_field(sp.x, sp.y, 0)
    
    analytical_fields={'Ez': resonant_cavity_ez_field}
    error_data = {
        'K': len(mesh_data['EToV']),
        'dt': driver.dt,
        'final_time': problem_data['n_steps'] * driver.dt,
        'time': [],
        'L2_error': {key: [] for key in analytical_fields}}
    t = 0.0

    for _ in range(int(problem_data['n_steps'])):
        error_data['time'].append(t)
        for field_name, analytical_fn in analytical_fields.items():
            uh = driver[field_name]
            ua = analytical_fn(sp.x, sp.y, t)
            l2_error = compute_L2_error(sp, uh, ua)
            error_data['L2_error'][field_name].append(l2_error)

        driver.step()
        t += driver.dt

    return error_data


def compute_rates(errors, epsilon=1e-14):
    """Remove erros pequenos, formata valores e calcula taxa média final."""
    row = []
    valid_errors = []

    for e in errors:
        if e < epsilon:
            row.append("—")
            valid_errors.append(np.nan)
        else:
            row.append(f"{e:.1E}")
            valid_errors.append(e)

    # Taxas de convergência
    rates = []
    for i in range(1, len(valid_errors)):
        if np.isnan(valid_errors[i-1]) or np.isnan(valid_errors[i]):
            rates.append(None)
        else:
            rates.append(np.log(valid_errors[i-1] / valid_errors[i]) / np.log(2))

    # Cálculo da taxa média dos últimos 2 valores
    avg_rate = (
        np.mean([r for r in rates[-2:] if r is not None])
        if any(r is not None for r in rates[-2:]) else None
    )

    row.append(f"{avg_rate:.1f}" if avg_rate else "—")
    return row, valid_errors


def plot_convergence_errors(problem_data, loglog_data, k_list, json_file, epsilon=1e-14):
    """Plota gráfico log-log dos erros L2 máximos."""

    # --- 1. Verifica se o JSON existe
    if not json_file.is_file():
        raise FileNotFoundError(f"Arquivo JSON não encontrado: {json_file}")

    # --- 3. Extrai pontos de referência
    fig_x, fig_y = extract_webdigitized_data(json_file)

    plt.figure(figsize=(8, 5))
    plt.scatter(fig_x, fig_y, s=8, marker='o', label='Hesthaven (2008)', color='black')
    for N, errors in loglog_data.items():
        ks = np.sqrt(np.array(k_list, dtype=float))
        errors = np.array(errors)
        mask = errors >= epsilon
        if np.any(mask):
            plt.loglog(ks[mask], errors[mask], marker='o', label=fr'Ordem $N={N}$')

    plt.xlabel(r'$K^{0.5}$')
    plt.ylabel('Erro $L^2$ Máximo')
    plt.title(fr"Fig. 6.9. Discrete $L^2$-error for $E_h^z$, obtained using a DG-FEM with {problem_data['flux_type']} flux.")
    plt.grid(True, which="both", ls="--", alpha=0.6)
    plt.xlim(1e0, 1e2)
    plt.ylim(1e-8, 1e0)
    plt.legend()
    plt.tight_layout()
    plt.show()


def run_convergence_rate_study(problem_data, json_name) -> pd.DataFrame:
    """
    Executa um estudo de convergência para o campo E.
    Retorna tabela de erros e gera gráfico.
    """
    epsilon = 1e-15
    N_list = [1, 2, 3, 4, 5]
    h_list = [2 / (2**n) for n in range(6)]

    table_data = {}
    loglog_data = {}

    print(f"\n🔎 Estudo de convergência | Fluxo: {problem_data['flux_type']}, N_steps: {problem_data['n_steps']}.")
    print("-------------------------------------------------------------------")

    for N in N_list:
        print(f"\n📐 Ordem polinomial N = {N}")
        raw_errors = []
        k_list = []
        for h in h_list:
            print(f"\n  ➤ Executando para K = {h} ...", end=" ")
            error_data = compute_error_data(problem_data, N, h)
            max_L2 = max(error_data['L2_error']['Ez'])
            raw_errors.append(max_L2)
            k_list.append(error_data['K'])
            print(f"Erro L2 máx = {max_L2:.2e}" if max_L2 >= epsilon else "Erro L2 máx ≈ 0")
            print(f"Run with {problem_data['n_steps']} steps and dt = {error_data['dt']} s. Final Time: {error_data['final_time']} s.")

        row, _ = compute_rates(raw_errors, epsilon)
        table_data[N] = row
        loglog_data[N] = raw_errors

    # Tabela
    col_labels = [str(k) for k in k_list] + ["Convergence rate"]
    df = pd.DataFrame.from_dict(table_data, orient='index', columns=col_labels)
    df.index.name = "N\\K"

    # Impressão
    print("\n✅ Estudo concluído. Tabela de erros L2 e taxas de convergência:")
    print(df)

    # Gráfico
    json_file = Path(__file__).parent.parent / 'examplesData' / 'inputs' / 'hesth_2D_cavity' / f"{json_name}.json"
    plot_convergence_errors(problem_data, loglog_data, k_list, json_file, epsilon=1e-14)
    return df


def main() -> None:
    """Função principal para execução do script."""
    clear_terminal()

    # PROBLEMA: Definições do problema de advecção linear 1D
    PROBLEM = {'name': 'hesth_2D_cavity_upwind',
        'folder_name': 'hesth_2D_cavity',
        'description': 'Teste de convergência do esquema DGTD bidimensional TMz. Hesthaven, p. 205',
        'json_name': 'hesthaven_fig_69a',
        'flux_type': 'Upwind',      # 'Upwind' or 'Centered'
        'cfl': 0.1,                 # Número de Courant-Friedrichs-Lewy
        'n_steps': 40,              # Tempo final da simulação
        'kx': 2*np.pi,              # Número de onda
        'L': 2,                     # Comprimento do domínio
        'n_order': 3                # Ordem de interpolação polinomial
    }
    
    # Test 0 - The solution
    single_test_solution(PROBLEM)    

    # Test 1 - Convergence study
    run_L2_error(PROBLEM)

    # Test 2 - Convergence rate study
    run_convergence_rate_study(PROBLEM, json_name=PROBLEM['json_name'])
    print("\n✅ Execução concluída com sucesso!")   


if __name__ == '__main__':
    main()
