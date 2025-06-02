##!/usr/bin/env python3
# -*- coding: utf-8 -*-
# pylint: disable=E0611, C0413, C0103, W0401, W0614, W0622
"""
══════════════════════════════════════════════════════════════════════════

ESTE SCRIPT UTILIZA CAMINHOS ABSOLUTOS E NÃO ALTERA O DIRETÓRIO DE TRABALHO

EXECUÇÃO:
cd C:\\git\\PyDG1D
conda activate pyDG1D
python examples\\cem_3p22.py

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


def main() -> None:
    """ Função principal para execução do script. """
    clear_terminal()

    sp = DG1D(
        n_order=3,
        mesh=Mesh1D(-1.0, 1.0, 20, boundary_label="SMA"),
        fluxType="Upwind"
    )

    driver = MaxwellDriver(sp, timeIntegratorType='LSERK4')

    # Senoidal initial condition
    # initialField = np.sin(2*np.pi*sp.x)

    # Gaussian initial condition
    x0 = 0.0  # Center of the Gaussian
    sigma = 0.1  # Standard deviation of the Gaussian
    initialField = np.exp(-((sp.x - x0) ** 2) / (2 * sigma ** 2))
    driver['E'][:] = initialField[:]
    #driver['H'][:] = initialField[:]

    # Cria a legenda com curvas não plotadas
    plt.figure()
    plt.plot([], [], 'b--', label=r'$E_y$ Field')
    plt.plot([], [], 'r-', label=r'$H_z$ Field')
    plt.legend(loc='upper right')  # Apenas uma vez, fora do loop

    t = 0.0
    for _ in range(100):
        plt.cla()  # Limpa os dados anteriores do gráfico

        # Atualiza os dados dos campos E e H
        plt.plot(sp.x, driver['E'], 'b--')
        plt.plot(sp.x, driver['H'], 'r-')

        plt.xlim(-1, 1)
        plt.ylim(-1.0, 1.0)
        plt.xlabel('x')
        plt.ylabel('Field Value')
        plt.title(f'Time: {t:.2f} seconds')
        plt.grid(True)

        # Adiciona a legenda novamente
        plt.plot([], [], 'b--', label=r'$E_y$ Field')
        plt.plot([], [], 'r-', label=r'$H_z$ Field')
        plt.legend(loc='upper right')
        plt.pause(0.05)
        driver.step()
        t += driver.dt
    plt.show()



if __name__ == '__main__':
    main()
