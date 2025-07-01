#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# pylint: disable=E0611, C0413, C0103
"""
cd C:\\Users\\adilt\\OneDrive\\05_GIT\\PYDG1D
conda activate pyDG1D
examples\jacobi_poly.py

------------------------------------------------------------
ESTE SCRIPT UTILIZA CAMINHOS ABSOLUTOS E NÃO ALTERA O DIRETÓRIO DE TRABALHO
------------------------------------------------------------

O script está estruturado para garantir que todos os arquivos e pastas usados
sejam acessados por caminhos absolutos, evitando problemas causados por mudanças
do diretório de trabalho (cwd).

PRINCIPAIS CONSIDERAÇÕES SOBRE DIRETÓRIOS:

1️⃣ Diretório raiz do projeto (CWD_ROOT):
    - CWD_ROOT = Path.cwd()
    - Definido automaticamente como o diretório atual quando o script é iniciado.
    - No seu ambiente, o diretório raiz é sempre:
      C:\\Users\\adilt\\OneDrive\\05_GIT\\PYDG1D
    - Todas as pastas e arquivos são referenciados em relação a esse diretório.

2️⃣ Estrutura esperada:
    - PYDG1D/
        ├── examples/
        │   └── jacobi_poly.py
        ├── maxwell/
        │   └── dg/
        │       ├── __init__.py
        │       ├── dg1d_tools.py
        │           └── jacobi_polynomial()
        │           └── jacobiGL()

Avaliação dos polinômios de Jacobi P_n^{(α,β)}(r) usando:
(1) Recorrência manual normalizada (jacobi_polynomial)
(2) Função do SciPy (eval_jacobi) com normalização

Compara os resultados e exibe o erro máximo para cada ordem.

Autor: Adilton PEREIRA
Data: 2025-05-19
"""

import os
import sys
import json
import numpy as np
from pathlib import Path
import matplotlib.pyplot as plt

from scipy.special import eval_jacobi
from scipy.special import gamma as gamma_func
from scipy.special import roots_jacobi

# Adiciona a raiz do projeto ao PYTHONPATH
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..\..')))

from maxwell.dg.dg1d_tools import jacobi_polynomial, jacobiGL, jacobi_gauss
from maxwell.utils import *

# Pastas principais
CASE_NAME   = 'jacobi_poly'
INPUTS      = (Path.cwd() / 'examplesData' / 'inputs' / CASE_NAME).resolve()

def extract_webdigitized_data(
    json_file_name: str,
    inputs_folder: Path = INPUTS
) -> tuple[np.ndarray, np.ndarray]:
    """
    Lê um arquivo JSON exportado pelo WebPlotDigitizer, extrai
    o primeiro dataset, ordena por X e retorna dois arrays: X e Y.
    
    Parâmetros
    ----------
    json_file_name : str
        Nome do arquivo JSON (deve estar em `inputs_folder`).
    inputs_folder : Path
        Pasta onde está o JSON.
    
    Retorno
    -------
    x, y : tuple de np.ndarray
        Vetores com as coordenadas ordenadas.
    """
    json_path = inputs_folder / json_file_name
    if not json_path.is_file():
        raise FileNotFoundError(f"JSON não encontrado: {json_path}")

    with open(json_path, 'r', encoding='utf-8') as f:
        obj = json.load(f)

    data_points = obj['datasetColl'][0]['data']
    # ordena pelo valor calibrado em X
    data_sorted = sorted(data_points, key=lambda pt: pt['value'][0])

    x = np.array([pt['value'][0] for pt in data_sorted])
    y = np.array([pt['value'][1] for pt in data_sorted])
    return x, y

def compute_gamma_n(
    n: int, alpha: float, beta: float
) -> float:
    """Calcula o fator de normalização γ_n para P_n^{(α,β)}(r)."""
    num = 2**(alpha + beta + 1) * gamma_func(n + alpha + 1) * gamma_func(n + beta + 1)
    den = (2 * n + alpha + beta + 1) * gamma_func(n + 1) * gamma_func(n + alpha + beta + 1)
    return num / den

def compare_jacobi_methods(
    alpha: float, beta: float, n_min: int = 1, n_max: int = 4
) -> None:
    """Compara os métodos de avaliação dos polinômios de Jacobi."""
    print("=== Comparação: Recorrência vs SciPy ===")
    for n_order in range(n_min, n_max + 1):

        # Método com recorrência
        r = jacobiGL(alpha, beta, n_order)
        Pn = jacobi_polynomial(r, alpha, beta, n_order)

        # Método SciPy com normalização
        gamma_n = compute_gamma_n(n_order, alpha, beta)
        Pn_scipy = eval_jacobi(n_order, alpha, beta, r) / np.sqrt(gamma_n)

        # Erro
        erro_max = np.max(np.abs(Pn - Pn_scipy))

        # Impressão formatada
        print(f"\nP_{n_order}^({alpha},{beta})(r):")
        print("r         =", np.round(r, 4))
        print("Manual    =", np.round(Pn, 6))
        print("SciPy     =", np.round(Pn_scipy, 6))
        print("gamma_n   =", round(gamma_n, 6))
        print(f"Erro máximo: {erro_max:.3e}")

def abramowitz_fig_22_1(
    json_file_name: str, alpha=1.5, beta=-0.5, n_min=1, n_max=5
) -> None:
    """
    Reproduz a Figura 22.1 do Abramowitz & Stegun usando a função eval_jacobi (não
    normalizados).
    A figura mostra os polinômios de Jacobi P_n^{(α,β)}(x) para diferentes ordens n.

    Handbook of Mathematical Functions de Abramowitz & Stegun, na:
    Seção: 22.2 — Orthogonality Relations
    """
    # 1) Extrai pontos de referência digitalizados
    abramo_x, abramo_y = extract_webdigitized_data(json_file_name)

    x = np.linspace(-1, 1, 400)
    plt.figure(figsize=(6, 8))

    # Plota os pontos digitalizados (referência)
    plt.scatter(abramo_x, abramo_y, s=3, marker='o', label='Abramowitz (1972)', color='black')

    for n in range(n_min, n_max + 1):
        #Pn = jacobi_polynomial(x, alpha, beta, n)
        Pn = eval_jacobi(n, alpha, beta, x)
        plt.plot(x, Pn, label=rf'$P_{{{n}}}^{{({alpha}, {beta})}}(x)$')

    plt.axhline(0, color='black', linewidth=0.5)
    plt.axvline(0, color='black', linewidth=0.5)
    plt.xlim(-1, 1)
    plt.ylim(-1.5, 3.5)
    plt.xlabel('$x$')
    plt.ylabel(r'$P_n^{(\alpha,\beta)}(x)$')
    plt.title(r'Jacobi Polynomials $P_n^{(\alpha,\beta)}(x)$'
               + f'\n$\\alpha = {alpha}, \\beta = {beta}, n = {n_min}(1){n_max}$')
    plt.grid(True, linestyle='--', alpha=0.3)
    plt.legend()
    plt.tight_layout()
    plt.show()

def compare_quadrature_nodes_weights(
    alpha: float, beta: float, N: int
) -> None:
    """
    Compara os pontos (x) e pesos (w) da quadratura de Gauss-Jacobi
    entre a implementação do PyDG1D (jacobi_gauss) e a do SciPy (roots_jacobi).

    Uma regra de quadratura de ordem 𝑁 integra exatamente todos os polinômios
    de grau até 2𝑁+1.

    Parâmetros:
    -----------
    alpha, beta : float
        Parâmetros dos polinômios de Jacobi.
    N : int
        Ordem da quadratura (polinômios de grau ≤ 2N+1 são integrados exatamente).
    """
    # PyDG1D (manual)
    x_pydg, w_pydg = jacobi_gauss(alpha, beta, N)

    # SciPy (nativa)
    x_scipy, w_scipy = roots_jacobi(N + 1, alpha, beta)

    # Ordena ambos os conjuntos (caso algum não esteja ordenado)
    x_pydg_sorted = x_pydg[np.argsort(x_pydg)]
    w_pydg_sorted = w_pydg[np.argsort(x_pydg)]
    x_scipy_sorted = x_scipy[np.argsort(x_scipy)]
    w_scipy_sorted = w_scipy[np.argsort(x_scipy)]

    print("\n=== Comparação: jacobi_gauss vs scipy.roots_jacobi ===")
    print(f"Ordem N = {N}  (grau polinômios integrados: 2N+1 = {2*N + 1})")
    print(f"alpha = {alpha}, beta = {beta}")
    print(f"{'x_pydg':>12} {'x_scipy':>12} {'Δx':>12}    |    {'w_pydg':>12} {'w_scipy':>12} {'Δw':>12}") # pylint: disable=C0301
    print("-" * 70)
    for i in range(N + 1):
        dx = abs(x_pydg_sorted[i] - x_scipy_sorted[i])
        dw = abs(w_pydg_sorted[i] - w_scipy_sorted[i])
        print(f"{x_pydg_sorted[i]:12.8f} {x_scipy_sorted[i]:12.8f} {dx:12.2e}    |    {w_pydg_sorted[i]:12.8f} {w_scipy_sorted[i]:12.8f} {dw:12.2e}") # pylint: disable=C0301

    # Erros máximos
    erro_max_x = np.max(np.abs(x_pydg_sorted - x_scipy_sorted))
    erro_max_w = np.max(np.abs(w_pydg_sorted - w_scipy_sorted))

    print(f"\nErro máximo nos nós: {erro_max_x:.3e}") 
    print(f"Erro máximo nos pesos: {erro_max_w:.3e}")

def compare_gauss_lobatto_nodes(
    alpha: float, beta: float, N: int
) -> None:
    """
    Compara os nós da quadratura Gauss-Lobatto-Jacobi obtidos por:
    (1) jacobiGL() do PyDG1D
    (2) Construção manual via extremos + raízes de P_{N-1}^{(α+1, β+1)}(x)

    Parâmetros:
    -----------
    alpha, beta : float
        Parâmetros dos polinômios de Jacobi.
    N : int
        Grau do polinômio interpolador (GL usa N+1 nós).
    """
    # 1. Pontos via PyDG1D
    r_pydg = jacobiGL(alpha, beta, N)

    # 2. Construção manual:
    # GL usa extremos + zeros do derivado → raízes de P_{N-1}^{(α+1, β+1)}
    if N == 1:
        r_manual = np.array([-1.0, 1.0])  # caso trivial
    else:
        r_interior, _ = roots_jacobi(N - 1, alpha + 1, beta + 1)
        r_manual = np.concatenate(([-1.0], r_interior, [1.0]))

    # 3. Ordenação para garantir alinhamento
    r_pydg_sorted = np.sort(r_pydg)
    r_manual_sorted = np.sort(r_manual)

    # 4. Comparação
    print("\n=== Comparação dos nós de Gauss-Lobatto ===")
    print(f"N = {N}  → {N+1} nós")
    print(f"alpha = {alpha}, beta = {beta}")
    print(f"{'r_pydg':>12} {'r_manual':>12} {'Δr':>12}")
    print("-" * 40)
    for i in range(N + 1):
        dr = abs(r_pydg_sorted[i] - r_manual_sorted[i])
        print(f"{r_pydg_sorted[i]:12.8f} {r_manual_sorted[i]:12.8f} {dr:12.2e}")

    # Erro máximo
    erro_max = np.max(np.abs(r_pydg_sorted - r_manual_sorted))
    print(f"\nErro máximo entre os nós: {erro_max:.3e}")

def main():
    """Função principal para execução do script."""
    clear_terminal()
    abramowitz_fig_22_1(json_file_name='abramowitz_fig_22_1.json')
    compare_jacobi_methods(alpha=0, beta=0)
    compare_quadrature_nodes_weights(alpha=0.0, beta=0.0, N=4)
    compare_gauss_lobatto_nodes(alpha=0.0, beta=0.0, N=4)

if __name__ == '__main__':
    main()
