#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# pylint: disable=E0611, C0413, C0103
""" 
VANDERMONDE_MATRICES.PY

Este script compara diferentes formas da matriz de Vandermonde associada aos
polinômios de Jacobi \( P_n^{(\alpha,\beta)}(r) \), conforme apresentado na
Figura 3.1 do livro de Hesthaven (2008), para grau de interpolação N = 6.

A matriz de Vandermonde é construída a partir da avaliação de funções base
em diferentes conjuntos de pontos nodais. Neste contexto, investigamos:

══════════════════════════════════════════════════════════════════════════
Figura 3.1 — Tipos de matrizes comparadas
══════════════════════════════════════════════════════════════════════════

┌────────┬─────────────────────┬────────────────────────────────────────    ┐
│ Figura │     Pontos ξ_i      │ Funções base ψ_n(r)                        │
├────────┼─────────────────────┼────────────────────────────────────────    ┤
│  (a)   │ Equidistantes       │ Monômios \( ψ_n = r^{n-1} \)               │
│  (b)   │ GLL (Gauss-Lobatto) │ Monômios \( ψ_n = r^{n-1} \)               │
│  (c)   │ Equidistantes       │ Jacobi ortonormalizados \( \tilde{P}_n \)  │
│  (d)   │ GLL (Gauss-Lobatto) │ Jacobi ortonormalizados \( \tilde{P}_n \)  │
└────────┴─────────────────────┴────────────────────────────────────────    ┘

Cada matriz é definida como:
    V_{ij} = ψ_j(r_i)

Para as formas ortonormalizadas, temos:
    ψ_j(r) = P_j^{(\alpha,\beta)}(r) / sqrt(γ_j)

em que:
    γ_j = fator de normalização conforme Apêndice A, Eq. (A.4)

══════════════════════════════════════════════════════════════════════════
FUNÇÕES DISPONÍVEIS NO SCRIPT:
══════════════════════════════════════════════════════════════════════════

1. compare_fig_31a(): 
   # (a) — Pontos equidistantes com base monomial.

2. compare_fig_31b(): 
   # (b) — Pontos GLL com base monomial.

3. compare_fig_31c(): 
   # (c) — Pontos equidistantes com base Jacobi ortonormalizada.

4. compare_fig_31d(): 
   # (d) — Pontos GLL com base Jacobi ortonormalizada.

Cada função compara uma matriz computada com a matriz correspondente do livro,
avaliando diferenças ponto a ponto (com arredondamento a 2 casas decimais).

══════════════════════════════════════════════════════════════════════════

ESTE SCRIPT UTILIZA CAMINHOS ABSOLUTOS E NÃO ALTERA O DIRETÓRIO DE TRABALHO

EXECUÇÃO:
    cd C:\\Users\\adilt\\OneDrive\\05_GIT\\PYDG1D
    conda activate pyDG1D
    python examples\\vandermonde_matrices.py

══════════════════════════════════════════════════════════════════════════

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
        │   └── vandermonde_matrices.py
        ├── maxwell/
        │   └── dg/
        │       ├── __init__.py
        │       ├── dg1d_tools.py
        │           └── jacobi_polynomial()
        │           └── jacobiGL()

Autor: Adilton Pereira
Data: 22/05/2025
"""

import os
import sys
import json
import numpy as np
from pathlib import Path
import matplotlib.pyplot as plt
from scipy.special import eval_jacobi
from scipy.special import eval_legendre
from scipy.special import gamma as gamma_func

# Adiciona a raiz do projeto ao PYTHONPATH
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..\..')))
from maxwell.dg.dg1d_tools import jacobiGL, vandermonde, differentiation_matrix

# Configurações do matplotlib
np.set_printoptions(precision=2, suppress=True)

# Pastas principais
CASE_NAME = 'vandermonde_matrices'
INPUTS = (Path.cwd() / 'examplesData' / 'inputs' / CASE_NAME).resolve()
ALPHA = 0
BETA  = 0

FIG31A = np.array([
    [1.00, -1.00,  1.00, -1.00,  1.00, -1.00,  1.00],
    [1.00, -0.67,  0.44, -0.30,  0.20, -0.13,  0.09],
    [1.00, -0.33,  0.11, -0.04,  0.01, -0.00,  0.00],
    [1.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00],
    [1.00,  0.33,  0.11,  0.04,  0.01,  0.00,  0.00],
    [1.00,  0.67,  0.44,  0.30,  0.20,  0.13,  0.09],
    [1.00,  1.00,  1.00,  1.00,  1.00,  1.00,  1.00]
])

FIG31B = np.array([
    [1.00, -1.00,  1.00, -1.00,  1.00, -1.00,  1.00],
    [1.00, -0.83,  0.69, -0.57,  0.48, -0.39,  0.33],
    [1.00, -0.47,  0.22, -0.10,  0.05, -0.02,  0.01],
    [1.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00],
    [1.00,  0.47,  0.22,  0.10,  0.05,  0.02,  0.01],
    [1.00,  0.83,  0.69,  0.57,  0.48,  0.39,  0.33],
    [1.00,  1.00,  1.00,  1.00,  1.00,  1.00,  1.00]
])

FIG31C = np.array([
    [0.71, -1.22,  1.58, -1.87,  2.12, -2.35,  2.55],
    [0.71, -0.82,  0.26,  0.49, -0.91,  0.72, -0.04],
    [0.71, -0.41, -0.53,  0.76,  0.03, -0.78,  0.49],
    [0.71,  0.00, -0.79,  0.00,  0.80,  0.00, -0.80],
    [0.71,  0.41, -0.53, -0.76,  0.03,  0.78,  0.49],
    [0.71,  0.82,  0.26, -0.49, -0.91, -0.72, -0.04],
    [0.71,  1.22,  1.58,  1.87,  2.12,  2.35,  2.55]
])

FIG31D = np.array([
    [0.71, -1.22,  1.58, -1.87,  2.12, -2.35,  2.55],
    [0.71, -1.02,  0.84, -0.35, -0.28,  0.81, -1.06],
    [0.71, -0.57, -0.27,  0.83, -0.50, -0.37,  0.85],
    [0.71,  0.00, -0.79,  0.00,  0.80,  0.00, -0.80],
    [0.71,  0.57, -0.27, -0.83, -0.50,  0.37,  0.85],
    [0.71,  1.02,  0.84,  0.35, -0.29, -0.81, -1.06],
    [0.71,  1.22,  1.58,  1.87,  2.12,  2.35,  2.55]
])

FIG34A = np.array([
    [-0.50,  0.50],
    [-0.50,  0.50]
])

FIG34B = np.array([
    [-1.50,  2.00, -0.50],
    [-0.50,  0.00,  0.50],
    [ 0.50, -2.00,  1.50]
])

FIG34C = np.array([
    [-5.00,  6.76, -2.67,  1.41, -0.50],
    [-1.24,  0.00,  1.75, -0.76,  0.26],
    [ 0.38, -1.34,  0.00,  1.34, -0.38],
    [-0.26,  0.76, -1.75,  0.00,  1.24],
    [ 0.50, -1.41,  2.67, -6.76,  5.00]
])

FIG34D = np.array([
    [-18.00,  24.35,  -9.75,   5.54,  -3.66,   2.59,  -1.87,   1.28,  -0.50],
    [ -4.09,   0.00,   5.79,  -2.70,   1.67,  -1.15,   0.82,  -0.56,   0.22],
    [  0.99,  -3.49,   0.00,   3.58,  -1.72,   1.08,  -0.74,   0.49,  -0.19],
    [ -0.44,   1.29,  -2.83,   0.00,   2.85,  -1.38,   0.86,  -0.55,   0.21],
    [  0.27,  -0.74,   1.27,  -2.66,   0.00,   2.66,  -1.27,   0.74,  -0.27],
    [ -0.21,   0.55,  -0.86,   1.38,  -2.85,   0.00,   2.83,  -1.29,   0.44],
    [  0.19,  -0.49,   0.74,  -1.08,   1.72,  -3.58,   0.00,   3.49,  -0.99],
    [ -0.22,   0.56,  -0.82,   1.15,  -1.67,   2.70,  -5.79,   0.00,   4.09],
    [  0.50,  -1.28,   1.87,  -2.59,   3.66,  -5.54,   9.75,  -24.35, 18.00]
])

def clear_terminal() -> None:
    """Limpa o terminal para melhor visualização."""
    os.system('cls' if os.name == 'nt' else 'clear')

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

def compare_fig_31a(N: int = 6) -> None:
    """Compara a matriz da Figura 3.1(a): base monomial nos pontos equidistantes."""
    print("\n=== Comparação: Figura 3.1(a) ===")

    # Pontos equidistantes no intervalo [-1, 1]
    r_eq = np.linspace(-1, 1, N + 1)

    # Matriz de Vandermonde com base monomial: V_{ij} = r_i^j
    vand_monomial = np.zeros((N + 1, N + 1))
    for n in range(N + 1):
        vand_monomial[:, n] = r_eq**n

    vand_rounded = np.round(vand_monomial, 2)
    diff = np.abs(vand_rounded - FIG31A)

    print(f"Erro máximo absoluto após arredondamento: {np.max(diff):.2e}")
    print("Matriz construída (monomial nos pontos equidistantes):")
    print(vand_rounded)

def compare_fig_31b(N: int = 6) -> None:
    """Compara a matriz da Figura 3.1(b): base monomial nos pontos GLL."""
    print("\n=== Comparação: Figura 3.1(b) ===")

    # Pontos de Gauss-Lobatto-Legendre (JacobiGL com α = β = 0)
    r_gll = jacobiGL(ALPHA, BETA, N)

    # Matriz de Vandermonde com base monomial: V_{ij} = r_i^j
    vand_monomial = np.zeros((N + 1, N + 1))
    for j in range(N + 1):
        vand_monomial[:, j] = r_gll**j

    vand_rounded = np.round(vand_monomial, 2)
    diff = np.abs(vand_rounded - FIG31B)

    print(f"Erro máximo absoluto após arredondamento: {np.max(diff):.2e}")
    print("Matriz construída (monomial nos pontos GLL):")
    print(vand_rounded)

def compare_fig_31c(N: int = 6) -> None:
    """Compara a matriz da Figura 3.1(c): base ortonormal nos pontos equidistantes."""
    print("\n=== Comparação: Figura 3.1(c) ===")

    r_eq = np.linspace(-1, 1, N + 1)  # Equidistantes
    vand_orthonormal = np.zeros((N + 1, N + 1))
    for n in range(N + 1):
        gamma_n = compute_gamma_n(n, 0, 0)  # fator de normalização
        Pn = eval_jacobi(n, 0, 0, r_eq)     # polinômio não normalizado
        vand_orthonormal[:, n] = Pn / np.sqrt(gamma_n)

    # Calcula diferença absoluta
    vand_rounded = np.round(vand_orthonormal, 2)
    diff = np.abs(vand_rounded - FIG31C)

    print(f"Erro máximo absoluto após arredondamento: {np.max(diff):.2e}")
    print("Matriz construída (ortonormal nos pontos equidistantes):")
    print(vand_rounded)

def compare_fig_31d(N: int = 6) -> None:
    """Compara a matriz da Figura 3.1(d): base ortonormal nos pontos LGL."""
    print("\n=== Comparação: Figura 3.1(d) ===")

    # Polinômios de Jacobi
    r = jacobiGL(ALPHA, BETA, N)
    vand_N6 = vandermonde(N, r)

    # Calcula diferença absoluta
    vand_rounded = np.round(vand_N6, 2)
    diff = np.abs(np.round(vand_N6, 2) - FIG31D)   
        
    print(f"Erro máximo absoluto após arredondamento: {np.max(diff):.2e}")
    print("Matriz construída (ortonormal nos pontos LGL):")
    print(vand_rounded)

def compare_fig_34a(N: int = 1) -> None:
    """Compara a matriz da Figura 3.4(a): Dr para N = 1."""
    print("\n=== Comparação: Figura 3.4(a) ===")

    # Polinômios de Jacobi
    r = jacobiGL(ALPHA, BETA, N)
    vandGrad_N1 = differentiation_matrix(N, r)

    # Calcula diferença absoluta
    vand_rounded = np.round(vandGrad_N1, 2)
    diff = np.abs(np.round(vandGrad_N1, 2) - FIG34A)   
        
    print(f"Erro máximo absoluto após arredondamento: {np.max(diff):.2e}")
    print("Matriz construída:")
    print(vand_rounded)

def compare_fig_34b(N: int = 2) -> None:
    """Compara a matriz da Figura 3.4(b): Dr para N = 2."""
    print("\n=== Comparação: Figura 3.4(b) ===")

    # Polinômios de Jacobi
    r = jacobiGL(ALPHA, BETA, N)
    vandGrad_N2 = differentiation_matrix(N, r)

    # Calcula diferença absoluta
    vand_rounded = np.round(vandGrad_N2, 2)
    diff = np.abs(np.round(vandGrad_N2, 2) - FIG34B)   
        
    print(f"Erro máximo absoluto após arredondamento: {np.max(diff):.2e}")
    print("Matriz construída:")
    print(vand_rounded)

def compare_fig_34c(N: int = 4) -> None:
    """Compara a matriz da Figura 3.4(c): Dr para N = 4."""
    print("\n=== Comparação: Figura 3.4(c) ===")

    # Polinômios de Jacobi
    r = jacobiGL(ALPHA, BETA, N)
    vandGrad_N4 = differentiation_matrix(N, r)

    # Calcula diferença absoluta
    vand_rounded = np.round(vandGrad_N4, 2)
    diff = np.abs(np.round(vandGrad_N4, 2) - FIG34C)   
        
    print(f"Erro máximo absoluto após arredondamento: {np.max(diff):.2e}")
    print("Matriz construída:")
    print(vand_rounded)

def compare_fig_34d(N: int = 8) -> None:
    """Compara a matriz da Figura 3.4(d): Dr para N = 8."""
    print("\n=== Comparação: Figura 3.4(d) ===")

    # Polinômios de Jacobi
    r = jacobiGL(ALPHA, BETA, N)
    vandGrad_N8 = differentiation_matrix(N, r)

    # Calcula diferença absoluta
    vand_rounded = np.round(vandGrad_N8, 2)
    diff = np.abs(np.round(vandGrad_N8, 2) - FIG34D)   
        
    print(f"Erro máximo absoluto após arredondamento: {np.max(diff):.2e}")
    print("Matriz construída:")
    print(vand_rounded)

def compute_determinants(N_values, node_type):
    """ Calcula o determinante da matriz de Vandermonde para diferentes N e tipos de nós."""        
    det_values = []
    for N in N_values:
        # Pontos de interpolação
        if node_type == 'equidistant':
            r = np.array([-1 + 2 * (i - 1) / N for i in range(1, N + 2)])
        elif node_type == 'gll':
            r = jacobiGL(ALPHA, BETA, N)
        else:
            raise ValueError("node_type deve ser 'equidistant' ou 'gll'.")

        # Monta a matriz de Vandermonde e calcula o determinante
        V = vandermonde(N, r)
        det_values.append(np.abs(np.linalg.det(V)))
    return det_values

def plot_fig_32(json_file_name: str, N_max: int = 34) -> None:
    """
    Exemplo 3.1 — Reproduz a Figura 3.2 do livro de Hesthaven,
    comparando o crescimento do determinante da matriz de Vandermonde
    construída com base ortonormal de polinômios de Legendre,
    para dois conjuntos de pontos:
        - Equidistant nodes (uniformemente espaçados)
        - LGL nodes (Legendre-Gauss-Lobatto)
    """
    # 1) Extrai pontos de referência digitalizados
    hest_x, hest_y = extract_webdigitized_data(json_file_name)

    # Geração dos dados
    N_vals = np.arange(2, N_max + 1)
    det_eq = compute_determinants(N_vals, node_type='equidistant')
    det_gll = compute_determinants(N_vals, node_type='gll')

    # Plotagem
    plt.figure(figsize=(8, 5))
    plt.scatter(hest_x, hest_y, s=3, marker='o', label='Hesthaven (2008)', color='black')
    plt.semilogy(N_vals, det_eq, '-.', label="Equidistant nodes")
    plt.semilogy(N_vals, det_gll, '--', label="LGL nodes")
    plt.xlabel("N", fontsize=12)
    plt.ylabel("|det(V)|", fontsize=12)
    plt.title("Fig. 3.2: Growth of the determinant of the Vandermonde matrix, V", fontsize=14)
    plt.xlim(2, N_max)
    xticks = list(range(2, 35, 4))  # [2, 6, 10, ..., 34]
    plt.xticks(xticks)
    plt.grid(True, which="both", linestyle="--", alpha=0.6)
    plt.legend()
    plt.tight_layout()
    plt.show()

def main():
    """Função principal para execução do script."""
    clear_terminal()
    print("Fig. 3.1. Entries of V for N = 6 and different choices of the basis, ψn(r), and evaluation points, ξi.")
    compare_fig_31a()
    compare_fig_31b()
    compare_fig_31c()
    compare_fig_31d()
    compare_fig_34a()
    compare_fig_34b()
    compare_fig_34c()
    compare_fig_34d()
    plot_fig_32(json_file_name='hesthaven_32.json')

if __name__ == '__main__':
    main()
