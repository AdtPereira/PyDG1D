#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# pylint: disable=E0611, C0413, C0103, W0401, W0614, W0622
"""
══════════════════════════════════════════════════════════════════════════

ESTE SCRIPT UTILIZA CAMINHOS ABSOLUTOS E NÃO ALTERA O DIRETÓRIO DE TRABALHO

EXECUÇÃO:
cd C:\\git\\PyDG1D
conda activate pyDG1D
python examples\\cem_3p21b.py

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


def monitor_local_energy(sp, driver, FinalTime: float, n_samples: int = 200) -> tuple[np.ndarray, np.ndarray]:
    """
    Monitoramento da energia local ao longo do tempo usando matriz de massa.
    """
    dt_sample = FinalTime / n_samples
    n_inner_steps = int(dt_sample / driver.dt)
    times = np.linspace(0, FinalTime, n_samples + 1)
    K = sp.mesh.number_of_elements()
    local_energy = np.zeros((K, n_samples + 1))

    for idx, t in enumerate(times):
        E = driver['E']
        H = driver['H']

        E_energy = sp.getEnergyPerElement(E)
        H_energy = sp.getEnergyPerElement(H)

        local_energy[:, idx] = 0.5 * (E_energy + H_energy)

        # Avança no tempo
        for _ in range(n_inner_steps):
            driver.step()

    return times, local_energy


def plot_local_energy(times, local_energy):
    """
    Gera um gráfico da energia local ao longo do tempo.

    Parâmetros
    ----------
    times : np.ndarray
        Vetor de tempos.
    local_energy : np.ndarray
        Energia local (shape: [K, len(times)]).
    """
    plt.figure(figsize=(10, 5))
    for k in range(local_energy.shape[0]):
        plt.plot(times, local_energy[k], label=f'Elemento {k+1}')
    plt.xlabel("Tempo [s]")
    plt.ylabel("Energia Local")
    plt.title("Evolução da Energia Local por Elemento")
    plt.grid(True, ls="--", alpha=0.6)
    plt.legend(ncol=2, fontsize=8)
    plt.tight_layout()
    plt.show()


def monitor_total_energy(sp, driver, FinalTime: float, n_samples: int = 200) -> tuple[np.ndarray, np.ndarray]:
    """
    Monitora a energia total armazenada no sistema ao longo do tempo.

    Parâmetros
    ----------
    sp : DG1D
        Objeto de discretização espacial com método getEnergy().
    driver : MaxwellDriver
        Driver da simulação que fornece os campos E e H.
    FinalTime : float
        Tempo final da simulação.
    n_samples : int
        Número de amostras no tempo.

    Retorno
    -------
    times : np.ndarray
        Vetor de instantes de tempo.
    total_energy : np.ndarray
        Vetor com a energia total no domínio para cada instante.
    """
    dt_sample = FinalTime / n_samples
    n_inner_steps = int(dt_sample / driver.dt)
    times = np.linspace(0.0, FinalTime, n_samples + 1)
    total_energy = np.zeros(n_samples + 1)

    for idx, t in enumerate(times):
        E = driver['E']
        H = driver['H']
        energy_E = sp.getEnergy(E)
        energy_H = sp.getEnergy(H)
        total_energy[idx] = 0.5 * (energy_E + energy_H)

        print(f"t = {t:.2f} s | energia total = {total_energy[idx]:.4e}")
        # Avança a simulação até a próxima amostragem
        for _ in range(n_inner_steps):
            driver.step()

    return times, total_energy


def plot_total_energy(times, total_energy):
    """
    Plota a energia total do sistema ao longo do tempo.
    """
    plt.figure(figsize=(8, 4))
    plt.plot(times, total_energy, 'k-', lw=2)
    plt.xlabel("Tempo [s]")
    plt.ylabel("Energia Total")
    plt.title("Evolução da Energia Total no Domínio")
    plt.ylim(0, 1)
    plt.grid(True, ls="--", alpha=0.6)
    plt.tight_layout()
    plt.show()


def main() -> None:
    """ Função principal para execução do script. """
    clear_terminal()

    # print("Iniciando estudo de convergência para o campo E...")
    # run_convergence_rate_study()

    # Monitoramento da energia local
    sp = DG1D(
        n_order=2,
        mesh=Mesh1D(0.0, 1.0, 16, boundary_label="PEC"),
        fluxType="Upwind")

    driver = MaxwellDriver(sp, timeIntegratorType='LSERK4')
    initialField = np.sin(kx*sp.x)
    driver['E'][:] = initialField[:]
    
    # times, local_energy = monitor_local_energy(
    #     sp,
    #     driver,
    #     FinalTime=20.0,
    #     n_samples=200
    # )

    # print(f"Monitoramento da energia local concluído. {len(times)} amostras coletadas.")
    # plot_local_energy(times, local_energy)

    time, total_energy = monitor_total_energy(sp, driver, FinalTime=1E3, n_samples=200)
    print(f"Monitoramento da energia total concluído. {len(time)} amostras coletadas.")
    plot_total_energy(time, total_energy)


if __name__ == '__main__':
    main()
