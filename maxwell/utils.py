""" Módulo de utilitários para o projeto Maxwell.
Este módulo contém funções auxiliares para manipulação de matrizes,
formatação de saídas e comparação com resultados do MATLAB.
"""

import os
import sys

# Adiciona a raiz do projeto ao PYTHONPATH
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.io import loadmat

# Configurações do matplotlib
np.set_printoptions(precision=4, suppress=True)


def clear_terminal() -> None:
    """Limpa o terminal para melhor visualização."""
    os.system('cls' if os.name == 'nt' else 'clear')


def format_matrix(matrix: np.ndarray, title: str = "", float_fmt: str = "%.4f") -> str:
    """
    Retorna uma string formatada da matriz com cabeçalho e valores bonitos.
    Pode ser usada para logs.

    Parâmetros
    ----------
    matrix : np.ndarray
        Matriz a ser formatada.
    title : str
        Título da matriz.
    float_fmt : str
        Formato dos valores.

    Retorno
    -------
    str
        String formatada.
    """
    df = pd.DataFrame(
        matrix, columns=[f"D{i+1}" for i in range(matrix.shape[1])])
    matrix_str = df.to_string(index=False, float_format=float_fmt)
    return f"\n{title}\n" + "-" * len(title) + f"\n{matrix_str}"


def compare_with_matlab(uh_py, mat_path, mat_var='u'):
    """ Compara a solução u da simulação Python com dados u do MATLAB. """
    # Caminho do arquivo .mat
    # mat_path = INPUTS / f"{PROBLEM['name']}.mat"
    if not mat_path.exists():
        raise FileNotFoundError(
            f'Arquivo MATLAB não encontrado: {mat_path}')

    # Carrega o arquivo .mat
    mat_data = loadmat(mat_path)
    if mat_var not in mat_data:
        raise KeyError(
            f"A variável '{mat_var} não foi encontrada em {mat_path}")

    # Remove dimensões extras, se houver
    mat_data = np.squeeze(mat_data[mat_var])

    # Verificação do tamanho
    if uh_py.shape != mat_data.shape:
        raise ValueError(
            f"Incompatibilidade: u_python = {uh_py.shape}, u_matlab = {mat_data.shape}")

    # Erro L2
    L2_norm = np.linalg.norm(uh_py - mat_data, ord=2)
    print(f'\nErro L2 entre soluções Python e MATLAB: {L2_norm:.3e}')

    return mat_data


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
    M = sp.mass            # (Np x Np)
    J = sp.jacobian        # (Np x K)
    K = sp.mesh.number_of_elements()
    errL2_local = np.zeros(K)

    for k in range(K):
        Jk = np.diag(J[:, k])         # (Np x Np)
        ek = err[:, k][:, np.newaxis] # (Np x 1)
        errL2_local[k] = (ek.T @ Jk @ M @ ek)[0, 0]
    return np.sqrt(np.sum(errL2_local))


def L2_error_E_field(sp, driver, analytical_E, n_steps):
    """
    Avalia a evolução da norma L2 do erro entre a solução numérica e analítica ao longo do tempo.

    Parâmetros
    ----------
    sp : DG1D
        Objeto do espaço DG contendo malha, pesos, etc.
    driver : MaxwellDriver
        Objeto com os campos numéricos e método .step().
    analytical_E : callable
        Função f(x, t) que retorna a solução analítica para o campo E_y.
    n_steps : int
        Número de passos de tempo a serem executados.

    Retorno
    -------
    dict
        Dicionário com chaves:
            'time'     → lista de instantes t
            'L2_error' → lista de ||u_h - u_a||_{L2}(t)
    """
    error_data = {'time': [], 'L2_error': []}
    t = 0.0

    for _ in range(n_steps):
        # Solução analítica reshape para shape (Np, K)
        ua = analytical_E(sp.x, t)
        uh = driver['E']
        l2_error = compute_L2_error(sp, uh, ua)
        error_data['time'].append(t)
        error_data['L2_error'].append(l2_error)
        driver.step()
        t += driver.dt
    return error_data


def L2_error_fields(sp, driver, analytical_fields: dict, n_steps: int) -> dict:
    """
    Avalia a evolução da norma L2 do erro entre as soluções numéricas e analíticas ao longo do tempo
    para múltiplos campos (E, H, etc.).

    Parâmetros
    ----------
    sp : DG1D
        Objeto do espaço DG contendo malha, pesos, etc.
    driver : MaxwellDriver
        Objeto com os campos numéricos e método .step().
    analytical_fields : dict
        Dicionário com chaves como 'E', 'H', etc., e valores sendo funções do tipo f(x, t).
    n_steps : int
        Número de passos de tempo a serem executados.

    Retorno
    -------
    dict
        Dicionário com:
            'time' → lista de instantes t
            'L2_error' → dicionário: nome do campo → lista de ||u_h - u_a||_{L2}(t)
    """
    error_data = {'time': [], 'L2_error': {key: [] for key in analytical_fields}}
    t = 0.0

    for _ in range(n_steps):
        error_data['time'].append(t)
        for field_name, analytical_fn in analytical_fields.items():
            uh = driver[field_name]
            ua = analytical_fn(sp.x, t)
            l2_error = compute_L2_error(sp, uh, ua)
            error_data['L2_error'][field_name].append(l2_error)

        driver.step()
        t += driver.dt

    return error_data
