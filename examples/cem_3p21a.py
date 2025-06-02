#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# pylint: disable=E0611, C0413, C0103, W0401, W0614, W0622
"""
══════════════════════════════════════════════════════════════════════════

ESTE SCRIPT UTILIZA CAMINHOS ABSOLUTOS E NÃO ALTERA O DIRETÓRIO DE TRABALHO

EXECUÇÃO:
cd C:\\git\\PyDG1D
conda activate pyDG1D
python examples\\cem_3p21a.py

══════════════════════════════════════════════════════════════════════════

O script está estruturado para garantir que todos os arquivos e pastas usados
sejam acessados por caminhos absolutos, evitando problemas causados por mudanças
do diretório de trabalho (cwd).

PRINCIPAIS CONSIDERAÇÕES SOBRE DIRETÓRIOS:

1️⃣ Diretório raiz do projeto (CWD_ROOT):
    - CWD_ROOT = Path.cwd()
    - Definido automaticamente como o diretório atual quando o script é iniciado.
    - No seu ambiente, o diretório raiz é sempre:
      C:\\git\\PYDG1D
    - Todas as pastas e arquivos são referenciados em relação a esse diretório.

2️⃣ Estrutura esperada:
    - PYDG1D/
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
Data: 30/05/2025
"""

import os
import sys

# Adiciona a raiz do projeto ao PYTHONPATH
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from maxwell.driver import *
from maxwell.dg.mesh1d import *
from maxwell.dg.dg1d import *
from maxwell.integrators.LSERK4 import *
from maxwell.utils import *

MODE = 2            # Modo de excitação da cavidade
c_ = 1              # Velocidade da luz no vácuo normalizada
eta_ = 1            # Impedância característica normalizada
kx = MODE * np.pi   # Número de onda


def analytical_solution_E(x, t):
    """ Solução analítica para o campo elétrico E. """
    return np.sin(kx * x) * np.cos(kx * c_ * t)


def analytical_solution_H(x, t):
    """ Solução analítica para o campo magnético H. """
    return - np.cos(kx * x) * np.sin(c_ * kx * t) / eta_


def max_L2_error(N, K, FinalTime, flux_type):
    """Executa a simulação para (N, K) e retorna o erro L2 máximo."""
    sp = DG1D(
        n_order=N,
        mesh=Mesh1D(0.0, 1.0, K, boundary_label="PEC"),
        fluxType=flux_type
    )

    driver = MaxwellDriver(sp, timeIntegratorType='LSERK4')
    driver['E'][:] = np.sin(kx * sp.x)

    error_data = L2_error_fields(
        sp,
        driver,
        analytical_fields={"E": analytical_solution_E},
        n_steps=int(FinalTime / driver.dt)
    )

    return max(error_data['L2_error']["E"])


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


def plot_convergence_errors(loglog_data, k_list, flux_type, FinalTime, epsilon=1e-14):
    """Plota gráfico log-log dos erros L2 máximos."""
    plt.figure(figsize=(8, 5))
    for N, errors in loglog_data.items():
        hs = 1 / np.array(k_list, dtype=float)
        errors = np.array(errors)
        mask = errors >= epsilon
        if np.any(mask):
            plt.loglog(hs[mask], errors[mask], marker='o', label=fr'Ordem $N={N}$')

    plt.xlabel(r'$h = 1/K$')
    plt.ylabel('Erro $L^2$ Máximo')
    plt.title('Estudo de Convergência (Erro L2 Máximo)\n' \
                f'Fluxo: {flux_type}, Tempo Final: {FinalTime:.0f} s')
    plt.grid(True, which="both", ls="--", alpha=0.6)
    plt.legend()
    plt.tight_layout()
    plt.show()


def run_convergence_rate_study(FinalTime=20.0, flux_type="Upwind") -> pd.DataFrame:
    """
    Executa um estudo de convergência para o campo E.
    Retorna tabela de erros e gera gráfico.
    """
    epsilon = 1e-15
    N_list = [1, 2, 4, 8]
    k_list = [2, 4, 8, 16, 32, 64]

    table_data = {}
    loglog_data = {}

    print(f"\n🔎 Estudo de convergência | Fluxo: {flux_type}, Tempo Final: {FinalTime:.2f} s")
    print("-------------------------------------------------------------------")

    for N in N_list:
        print(f"\n📐 Ordem polinomial N = {N}")
        raw_errors = []
        for K in k_list:
            print(f"  ➤ Executando para K = {K} ...", end=" ")
            max_L2 = max_L2_error(N, K, FinalTime, flux_type)
            raw_errors.append(max_L2)
            print(f"Erro L2 máx = {max_L2:.2e}" if max_L2 >= epsilon else "Erro L2 máx ≈ 0")

        row, _ = compute_rates(raw_errors, epsilon)
        table_data[N] = row
        loglog_data[N] = raw_errors

    # Tabela
    col_labels = [str(k) for k in k_list] + ["Convergence rate"]
    df = pd.DataFrame.from_dict(table_data, orient='index', columns=col_labels)
    df.index.name = "N\\K"

    # Impressão
    print("\n✅ Estudo concluído. Tabela de erros L2 e taxas de convergência:")
    print(f"Fluxo: {flux_type}, Tempo Final: {FinalTime:.2f} s")
    print(df)

    # Gráfico
    plot_convergence_errors(loglog_data, k_list, flux_type, FinalTime, epsilon=1e-14)
    return df


def main() -> None:
    """ Função principal para execução do script. """
    clear_terminal()

    print("Iniciando estudo de convergência para o campo E...")
    run_convergence_rate_study()


if __name__ == '__main__':
    main()
