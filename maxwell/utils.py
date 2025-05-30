""" Módulo de utilitários para o projeto Maxwell.
Este módulo contém funções auxiliares para manipulação de matrizes,
formatação de saídas e comparação com resultados do MATLAB.
"""

import os
import numpy as np
import pandas as pd
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


def compute_L2_error(sp, u_h, ua):
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
    err = ua - u_h
    M = sp.mass            # (Np x Np)
    J = sp.jacobian        # (Np x K)
    K = sp.mesh.number_of_elements()
    errL2_local = np.zeros(K)

    for k in range(K):
        Jk = np.diag(J[:, k])         # (Np x Np)
        ek = err[:, k][:, np.newaxis] # (Np x 1)
        errL2_local[k] = (ek.T @ Jk @ M @ ek)[0, 0]

    return np.sqrt(np.sum(errL2_local))


