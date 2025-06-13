#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# pylint: disable=E0611, C0413, C0103, W0401, W0614, W0622
"""
══════════════════════════════════════════════════════════════════════════

ESTE SCRIPT UTILIZA CAMINHOS ABSOLUTOS E NÃO ALTERA O DIRETÓRIO DE TRABALHO

EXECUÇÃO:
cd C:\\git\\PyDG1D
conda activate pyDG1D
python examples\\cem_4p8.py

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

# Adiciona a raiz do projeto ao PYTHONPATH
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..\..')))

import numpy as np
import pandas as pd
from pathlib import Path
import matplotlib.pyplot as plt
from matplotlib.tri import Triangulation

from maxwell.dg.dg2d import *
from maxwell.driver import *
from maxwell.integrators.LSERK4 import *
from mesher.create_mesh import *
from maxwell.utils import *


BOUNDARY = [{'tag': 101, 'type': 'Dirichlet', 'value': 0.0, 'name': 'Ez_0'}]
MATERIAL = [{'tag': 201, 'name': 'free_space', 'relative_magnetic_permeability': 1, 'relative_electric_permittivity': 1}]   
INFO_GRAPH = {'cell': False, 'nodes': False, 'edges': False, 'edges_numb': False, 'filepath': 'examplesData/inputs/cem_4p8/cem_4p8.svg'}
PROBLEM = {'name': 'cem_4p8', 'folder_name': 'cem_4p8', 'description': 'Teste de convergência do esquema DGTD bidimensional TMz.',
    'bc': "PEC",                # Condição de contorno: 'PEC'  or 'Periodic
    'flux_type': 'Centered',      # 'Upwind' or 'Centered'
    'cfl': 1E-1,                  # Número de Courant-Friedrichs-Lewy
    'm': 1,                     # Número de modo
    'n': 1,                     # Número de modo
    'L': np.pi,                 # Dimensão total do domínio
    'n_order': 3                # Ordem de interpolação polinomial
}


class SpectralAnalyzer:
    def __init__(self, PROBLEM):
        self.m = PROBLEM['m']
        self.n = PROBLEM['n']
        self.L = PROBLEM['L']
        kx = self.m * np.pi / self.L
        ky = self.n * np.pi / self.L
        self.w = np.sqrt(kx**2 + ky**2)
        self.t_final = 3 * (2 * np.pi / self.w)

    def resonant_cavity_ez_field(self, x, y, t):
        return np.cos(x) * np.cos(y) * np.cos(self.w*t)


    def resonant_cavity_hx_field(self, x, y, t):
        return np.cos(x) * np.sin(y) * np.sin(self.w*t) / self.w


    def resonant_cavity_hy_field(self, x, y, t):
        return - np.sin(x) * np.cos(y) * np.sin(self.w*t) / self.w


    def plot_analytical_field(self, Nx=200, Ny=200):
        """
        Plota a solução analítica Ez sobre o domínio físico [-π, π]^2
        em t = 3T, com dois subplots: mapa de cores 2D e superfície 3D.
        """
        # Domínio
        L = self.L
        T = 2 * np.pi / self.w
        x = np.linspace(-L/2, L/2, Nx)
        y = np.linspace(-L/2, L/2, Ny)
        X, Y = np.meshgrid(x, y)

        # Tempo final = 3 períodos
        t_final = 3 * T

        # Avaliação da solução
        Ez = self.resonant_cavity_ez_field(X, Y, t_final)

        # Criação dos subplots
        fig = plt.figure(figsize=(12, 5))
        fig.suptitle(rf'Solução Analítica $E_z(x, y, t = 3T = {t_final:.2f} s)$', fontsize=14)

        # Subplot 1: Heatmap 2D
        ax1 = fig.add_subplot(1, 2, 1)
        pcm = ax1.pcolormesh(X, Y, Ez, shading='auto', cmap='viridis')
        ax1.set_xlabel('x')
        ax1.set_ylabel('y')
        ax1.axis('equal')
        fig.colorbar(pcm, ax=ax1, label=r'$E_z$')

        # Subplot 2: Superfície 3D
        ax2 = fig.add_subplot(1, 2, 2, projection='3d')
        surf = ax2.plot_surface(X, Y, Ez, cmap='viridis', edgecolor='none')
        ax2.set_xlabel('x')
        ax2.set_ylabel('y')
        ax2.set_zlabel(r'$E_z$')
        fig.colorbar(surf, ax=ax2, shrink=0.6, aspect=12)

        plt.tight_layout()
        plt.show()


def single_test_solution(PROBLEM) -> None:
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
    sa = SpectralAnalyzer(PROBLEM)
    # Criar a malha retangular com Gmsh
    mesh_data = mesh_rectangular_domain(PROBLEM, BOUNDARY, MATERIAL, h=1, view_mesh=False, mesh_info=False) 

    # Visualizar a malha criada
    # plot_triangular_mesh(INFO_GRAPH, mesh_data)

    # Solução analítica
    # sa.plot_analytical_field()

    # Definir a discretização espacial usando DG2D
    sp = Maxwell2D(
        n_order=PROBLEM['n_order'],
        mesh=Mesh2D(
            vx=mesh_data['VX'],
            vy=mesh_data['VY'],
            EToV=mesh_data['EToV'],
            boundary_label=PROBLEM['bc']),
        fluxType=PROBLEM['flux_type'])
    
    # Initialize the solver
    driver = MaxwellDriver(sp, CFL=PROBLEM['cfl'])
    driver['Ez'][:] = sa.resonant_cavity_ez_field(sp.x, sp.y, 0)
    
    # Set the time integrator
    driver.run_until(sa.t_final)
    n_steps = int(np.ceil(sa.t_final / driver.dt))
    print(f"\n🌐 Passo de tempo da simulação: {driver.dt:.3f} s")
    print(f"🌐 Número de passos de tempo: {n_steps}")
    print(f"🌐 Tempo final da simulação: {driver.timeIntegrator.time:.3f} s")
    
    ez_expected = sa.resonant_cavity_ez_field(sp.x, sp.y, driver.timeIntegrator.time)
    R = np.corrcoef(ez_expected, driver['Ez'])
    print(f"🌐 Correlação entre solução esperada e solução numérica: {R[0,1]:.2f}")

    # Solução numérica e analítica finais
    # print(f"\n🌐 Solução numérica final Ez:\n{driver['Ez']}\n")
    # print(f"🌐 Solução analítica final Ez:\n{ez_expected}\n")

    # Plota a solução numérica Ez diretamente nos pontos modais
    plot_global_solution(sp, driver['Ez'], field_label='E_z', show_mesh=True)


def plot_global_solution(sp, field, field_label='u', show_mesh=True):
    """
    Plota a solução DG diretamente nos pontos modais do método, sem interpolação.

    Parâmetros
    ----------
    sp : Maxwell2D
        Objeto com malha e campos do método DG.
    field : ndarray (Np x K)
        Solução numérica DG (Ez, Hx ou Hy).
    field_label : str
        Nome do campo a ser exibido no gráfico.
    show_mesh : bool
        Se True, sobrepõe as bordas dos elementos triangulares.
    """
    # Reorganiza os dados em arrays planos
    x = sp.x.ravel(order='F')  # tamanho Np*K
    y = sp.y.ravel(order='F')
    uh = field.ravel(order='F')

    r, s = xy_to_rs(*set_nodes_in_equilateral_triangle(sp.n_order))
    tri_local = Triangulation(r, s).triangles

    # Constrói triangulação global
    Np, K = field.shape
    triangles_all = []
    for k in range(K):
        offset = k * Np
        tri_k = tri_local + offset
        triangles_all.extend(tri_k)
    tri_global = Triangulation(x, y, np.array(triangles_all))

    # --- Gráficos ---
    fig = plt.figure(figsize=(14, 6))
    ax1 = fig.add_subplot(1, 2, 1)
    ax2 = fig.add_subplot(1, 2, 2, projection='3d')

    # Contorno 2D
    tcf = ax1.tricontourf(tri_global, uh, levels=100, cmap='viridis')
    fig.colorbar(tcf, ax=ax1, label=fr'{field_label}(x,y)')
    ax1.set_title(fr'Solução DG ${field_label}(x, y)$ nos pontos modais')
    ax1.set_xlabel("x")
    ax1.set_ylabel("y")
    ax1.set_aspect("equal")

    if show_mesh:
        VX, VY = sp.mesh.vx, sp.mesh.vy
        EToV = sp.mesh.EToV
        for tri in EToV:
            px = [VX[tri[0]], VX[tri[1]], VX[tri[2]], VX[tri[0]]]
            py = [VY[tri[0]], VY[tri[1]], VY[tri[2]], VY[tri[0]]]
            ax1.plot(px, py, color='k', linewidth=0.5)

    # Superfície 3D
    ax2.plot_trisurf(tri_global, uh, cmap='viridis', linewidth=0.2)
    ax2.set_title(fr'Solução DG ${field_label}(x, y)$ - Superfície')
    ax2.set_xlabel("x")
    ax2.set_ylabel("y")
    ax2.set_zlabel(fr'{field_label}(x, y)')
    ax2.view_init(elev=35, azim=-120)
    plt.tight_layout()
    plt.show()


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


def run_L2_error(PROBLEM):
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
    sa = SpectralAnalyzer(PROBLEM)
    t_final = sa.t_final
    h_list = [2]
    for h in h_list:
        # Atualiza o número de elementos na malha
        mesh_data = mesh_rectangular_domain(PROBLEM, BOUNDARY, MATERIAL, h, auto_save=False)

        # Generate Spatial Discretization
        sp = Maxwell2D(
            n_order=PROBLEM['n_order'],
            mesh=Mesh2D(
                vx=mesh_data['VX'],
                vy=mesh_data['VY'],
                EToV=mesh_data['EToV'],
                boundary_label="PEC"),
            fluxType=PROBLEM['flux_type'])
        
        # Initialize the solver
        driver = MaxwellDriver(sp, CFL=PROBLEM['cfl'])
        driver['Ez'][:] = sa.resonant_cavity_ez_field(sp.x, sp.y, 0)
        
        analytical_fields={'Ez': sa.resonant_cavity_ez_field}
        error_data = {'time': [], 'L2_error': {key: [] for key in analytical_fields}}
        t = 0.0

        for _ in range(int(t_final / driver.dt)):
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


def compute_error_data(PROBLEM, N, h):
    """Executa a simulação para (N, K) e retorna o erro L2 máximo."""
    sa = SpectralAnalyzer(PROBLEM)
    # Atualiza o número de elementos na malha
    mesh_data = mesh_rectangular_domain(PROBLEM, BOUNDARY, MATERIAL, h, auto_save=False)

    # Generate Spatial Discretization
    sp = Maxwell2D(
        n_order=N,
        mesh=Mesh2D(
            vx=mesh_data['VX'],
            vy=mesh_data['VY'],
            EToV=mesh_data['EToV'],
            boundary_label=PROBLEM['bc']),
        fluxType=PROBLEM['flux_type'])

    # Initialize the solver
    driver = MaxwellDriver(sp, CFL=PROBLEM['cfl'])
    driver['Ez'][:] = sa.resonant_cavity_ez_field(sp.x, sp.y, 0)
    n_steps = int(np.ceil(sa.t_final / driver.dt))

    analytical_fields={'Ez': sa.resonant_cavity_ez_field}
    error_data = {
        'K': len(mesh_data['EToV']),
        'dt': driver.dt,
        'n_steps': n_steps,
        'time': [],
        'L2_error': {key: [] for key in analytical_fields}}

    t = 0.0
    for _ in range(n_steps):
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


def plot_convergence_errors(PROBLEM, loglog_data, k_list, epsilon=1e-14):
    """Plota gráfico log-log dos erros L2 máximos."""

    plt.figure(figsize=(8, 5))
    for N, errors in loglog_data.items():
        ks = np.sqrt(np.array(k_list, dtype=float))
        errors = np.array(errors)
        mask = errors >= epsilon
        if np.any(mask):
            plt.loglog(ks[mask], errors[mask], marker='o', label=fr'Ordem $N={N}$')

    plt.xlabel(r'$K^{0.5}$')
    plt.ylabel('Erro $L^2$ Máximo')
    plt.title(fr"DG-FEM Discrete $L^2$-error for $E_h^z$, for {PROBLEM['bc']} and {PROBLEM['flux_type']} flux.")
    plt.grid(True, which="both", ls="--", alpha=0.6)
    plt.xlim(1e0, 1e2)
    plt.ylim(1e-10, 1e0)
    plt.legend()
    plt.tight_layout()
    plt.show()


def run_convergence_rate_study(PROBLEM) -> pd.DataFrame:
    """
    Executa um estudo de convergência para o campo E.
    Retorna tabela de erros e gera gráfico.
    """
    sa = SpectralAnalyzer(PROBLEM)
    epsilon = 1e-15
    N_list = [1, 2, 3, 4, 5]
    h_list = [4 / (2**n) for n in range(5)] # [4, 2, 1, 0.5, 0.25]

    table_data = {}
    loglog_data = {}

    print(f"\n🔎 Estudo de convergência | Fluxo: {PROBLEM['flux_type']}, t_final: {sa.t_final:.2f}.")
    print("-------------------------------------------------------------------")

    for N in N_list:
        print(f"\n📐 Ordem polinomial N = {N}")
        raw_errors = []
        k_list = []
        for h in h_list:
            print(f"\n  ➤ Executando para K = {h} ...", end=" ")
            error_data = compute_error_data(PROBLEM, N, h)
            max_L2 = max(error_data['L2_error']['Ez'])
            raw_errors.append(max_L2)
            k_list.append(error_data['K'])
            print(f"Erro L2 máx = {max_L2:.2e}" if max_L2 >= epsilon else "Erro L2 máx ≈ 0")
            print(f"Run with {error_data['n_steps']} steps and dt = {error_data['dt']} s.")

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
    plot_convergence_errors(PROBLEM, loglog_data, k_list, epsilon=1e-14)
    return df


def main() -> None:
    """Função principal para execução do script."""
    clear_terminal()

    # Test 0 - The solution
    single_test_solution(PROBLEM)    

    # Test 1 - Convergence study
    # run_L2_error(PROBLEM)

    # Test 2 - Convergence rate study
    # run_convergence_rate_study(PROBLEM)
    # print("\n✅ Execução concluída com sucesso!")


if __name__ == '__main__':
    main()
