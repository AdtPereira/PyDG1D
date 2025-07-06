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
import json
from typing import List, Optional
from scipy.io import loadmat

# Configurações do matplotlib
np.set_printoptions(precision=4, suppress=True)


def clear_terminal() -> None:
    """Limpa o terminal para melhor visualização."""
    os.system('cls' if os.name == 'nt' else 'clear')


def print_3d_matrices(Matriz3D: np.ndarray, elements: Optional[List[int]] = None, title: str = "") -> None:
    """
    Imprime fatias de vmapM_3D usando DataFrames do pandas para uma exibição elegante.

    Permite a seleção de quais elementos (fatias ao longo do último eixo) imprimir.

    Argumentos:
        vmapM_3D (np.ndarray): A matriz com shape (Nfp, Nfaces, K).
        title (str): O título a ser exibido no cabeçalho.
        elementos_a_plotar (Optional[List[int]], opcional): 
            Uma lista de índices (k) dos elementos a serem plotados.
            Se for None, todos os elementos (de 0 a K-1) serão plotados. 
            O padrão é None.
    """
    if Matriz3D.ndim != 3:
        print("Erro: A matriz não é 3D.")
        return
        
    Nfp, Nfaces, K = Matriz3D.shape

    # Determina quais elementos imprimir
    if elements is None:
        idx_to_print = range(K)  # Imprime todos os elementos
    elif all(isinstance(el, int) for el in elements):
        idx_to_print = elements
    else:
        print("Erro: 'elements' deve ser uma lista de índices ou None.")
        return

    print("\n" + "="*40)
    print(f"       {title} ")
    print("="*40)

    printed_elements = 0
    for k in idx_to_print:
        # Verificação de segurança: garante que o índice está dentro dos limites
        if not 0 <= k < K:
            print(f"\nAviso: O índice de elemento {k} está fora do intervalo [0, {K-1}] e será ignorado.")
            continue

        print(f"\n--- Elemento {k} ---")

        slice_k = Matriz3D[:, :, k]

        # Cria um DataFrame do pandas a partir da fatia 2D
        df = pd.DataFrame(slice_k)

        # Adiciona rótulos para as colunas e o índice (linhas)
        df.columns = [f'Face {i}' for i in range(Nfaces)]
        df.index.name = "Face Node"

        print(df)
        printed_elements += 1

    if printed_elements == 0 and elements is not None:
        print("\nNenhum elemento válido foi encontrado para impressão.")


def format_matrix(
    matrix: np.ndarray, 
    title: str = "", 
    elements: Optional[List[int]] = None,
    float_fmt: str = "%.4f"
) -> str:
    """
    Retorna uma string formatada da matriz com cabeçalho e valores.
    Permite a seleção de quais linhas da matriz formatar.

    Parâmetros
    ----------
    matrix : np.ndarray
        Matriz a ser formatada.
    title : str
        Título da matriz.
    float_fmt : str
        Formato dos valores de ponto flutuante.
    linhas_a_formatar : Optional[List[int]], opcional
        Uma lista de índices das linhas a serem incluídas na formatação.
        Se for None, todas as linhas da matriz serão usadas. O padrão é None.

    Retorno
    -------
    str
        String formatada.
    """
    if matrix.ndim != 2:
        return f"\n{title}\n" + "-" * len(title) + "\nErro: A matriz de entrada não é 2D."

    num_rows = matrix.shape[0]
    matriz_para_df = matrix

    if elements is not None:
        # Filtra os índices para garantir que estão dentro dos limites da matriz
        indices_validos = [i for i in elements if 0 <= i < num_rows]
        
        # Se após a filtragem não sobrar nenhum índice válido, retorna uma mensagem
        if not indices_validos:
            return f"\n{title}\n" + "-" * len(title) + "\nNenhuma linha válida foi selecionada para formatação."
            
        # Seleciona apenas as linhas válidas da matriz original usando a indexação do NumPy
        matriz_para_df = matrix[indices_validos, :]

    # Cria o DataFrame a partir da matriz (completa ou fatiada)
    if title == "EToV":
        columns = [f"Vertex {i}" for i in range(matrix.shape[1])]
    elif title == "EToE" or title == "EToF":
        columns = [f"Face {i}" for i in range(matrix.shape[1])]
    else:
        # Para outros títulos, usa D1, D2, ..., Dn
        columns = [f"D{i}" for i in range(matrix.shape[1])]

        
    df = pd.DataFrame(
        matriz_para_df, columns=columns, index=[f"t{i}" for i in indices_validos]
    )
    
    # Formata o DataFrame como string
    matrix_str = df.to_string(index=True, float_format=float_fmt) # Mudei para index=True para ver os índices originais
    
    # Monta a string final com o título
    header = f"\n{title} " + f"(dim: {matrix.shape})\n" + "-" * 3*len(title)

    # Imprimir resultado
    print(f"{header}\n{matrix_str}")

    return f"{header}\n{matrix_str}"


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


def extract_webdigitized_data(json_path: str) -> tuple[np.ndarray, np.ndarray]:
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
