#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# pylint: disable=E0611, C0413, C0103, W0401, W0614, W0622
"""
══════════════════════════════════════════════════════════════════════════

ESTE SCRIPT UTILIZA CAMINHOS ABSOLUTOS E NÃO ALTERA O DIRETÓRIO DE TRABALHO

EXECUÇÃO:
cd C:\\git\\PyDG1D
conda activate pyDG1D
python examples\\cem_3p21.py

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


def animate_fields(sp, driver, analytical_E=None, analytical_H=None, n_frames=150):
    """
    Anima a evolução dos campos E_y e H_z, com comparações analíticas opcionais.

    Parâmetros
    ----------
    sp : DG1D
        Espaço de discretização.
    driver : MaxwellDriver
        Driver contendo os campos numéricos e método .step().
    analytical_E : callable or None
        Função f(x, t) para solução analítica de E_y. Se None, não plota.
    analytical_H : callable or None
        Função f(x, t) para solução analítica de H_z. Se None, não plota.
    n_frames : int
        Número de quadros da animação.
    """

    t = 0.0
    for _ in range(n_frames):
        plt.cla()  # Limpa os dados anteriores do gráfico

        # Atualiza os dados numéricos
        plt.plot(sp.x, driver['E'], 'b--')
        plt.plot(sp.x, driver['H'], 'r-')

        # Atualiza os dados analíticos, se fornecidos
        if analytical_E is not None:
            E_analytical = analytical_E(sp.x, t)
            plt.plot(sp.x, E_analytical, 'k.')
        if analytical_H is not None:
            H_analytical = analytical_H(sp.x, t)
            plt.plot(sp.x, H_analytical, 'k.')

        plt.xlim(0, 1)
        plt.ylim(-1.5, 1.5)
        plt.xlabel('x')
        plt.ylabel('Field Value')
        plt.title(f'Time: {t:.2f} seconds')
        plt.grid(True)

        # Adiciona a legenda novamente
        plt.plot([], [], 'b--', label=r'$E_y$ Field')
        plt.plot([], [], 'r-', label=r'$H_z$ Field')
        plt.plot([], [], 'k.', label=r'$Analytical')
        plt.legend(loc='upper right')
        plt.pause(0.05)
        driver.step()
        t += driver.dt


def plot_L2_error_evolution(error_data):
    """
    Plota a evolução temporal da norma L2 do erro.

    Parâmetros
    ----------
    error_data : dict
        Dicionário com as chaves:
            'time'     → lista de tempos t
            'L2_error' → lista de ||u_h - u_a||_{L2}(t)
    """
    time = error_data['time']
    error = error_data['L2_error']

    plt.figure(figsize=(8, 4))
    plt.plot(
        time, error, 'k-o', linewidth=1.5, markersize=4, label=r'$\|E_h - E_a\|_{L^2}$')
    plt.xlabel('Tempo')
    plt.ylabel('Erro $L^2$')
    plt.title('Evolução Temporal do Erro $L^2$')
    plt.grid(True, linestyle='--', alpha=0.7)
    plt.legend()
    plt.tight_layout()
    plt.show()


def main() -> None:
    """ Função principal para execução do script. """
    clear_terminal()

    sp = DG1D(
        n_order=3,
        mesh=Mesh1D(0.0, 1.0, 16, boundary_label="PEC"),
        fluxType="Centered"
    )

    driver = MaxwellDriver(sp, timeIntegratorType='LSERK4')
    initialField = np.sin(kx*sp.x)
    driver['E'][:] = initialField[:]

    # Animação dos campos E e H com soluções analíticas
    animate_fields(
        sp, driver,
        analytical_E=analytical_solution_E,
        analytical_H=analytical_solution_H,
        n_frames=150)

    # Cálculo do erro L2 entre a solução numérica e analítica
    # error_E_field = L2_error_E_field(
    #     sp, driver,
    #     analytical_E=analytical_solution_E,
    #     n_steps=150)

    # Valor máximo do Erro L2
    # print(f"Erro L2 máximo para o campo Elétrico: {max(error_E_field['L2_error']):.4e}")

    # Plotando o erro L2 ao longo do tempo
    # plot_L2_error_evolution(error_E_field)

if __name__ == '__main__':
    main()
