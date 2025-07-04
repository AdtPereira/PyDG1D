#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# pylint: disable=E0611, C0413, C0103, W0401, W0614, W0622
"""
══════════════════════════════════════════════════════════════════════════

ESTE SCRIPT UTILIZA CAMINHOS ABSOLUTOS E NÃO ALTERA O DIRETÓRIO DE TRABALHO

EXECUÇÃO:
cd C:\\git\\PyDG1D
conda activate pyDG1D
python examples\\cem_5\\upml.py

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
import matplotlib.pyplot as plt

# Adiciona a raiz do projeto ao PYTHONPATH
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..\..')))

from maxwell.dg.dg3d import *
from maxwell.driver import *
from maxwell.utils import *
from maxwell.integrators.LSERK4 import *
from maxwell.spectralAnalyzer import UniaxialPml

from mesher.create_mesh import *
from mesher.plot_mesh import *


PROBLEM = {
    'DIM': 3,
    'description': 'Teste da UPML do esquema DGTD tridimensional.',
    'name': 'upml',
    'folder': 'cem_5',
    'dg': {
        'n_order': 4,               # Ordem de interpolação polinomial
        'flux_type': 'Upwind',      # 'Upwind' or 'Centered'
        'cfl': 1.0,                 # Número de Courant-Friedrichs-Lewy
        'bc': "SMA",                # Condição de contorno: 'PEC', 'SMA'  or 'Periodic
        't_final': 10.0,            # Tempo final da simulação
    },
    'domain': {
        'type': 'cubic',            # Tipo de domínio: 'rectangle' ou 'cubic'
        'h': 1.0,                   # Tamanho máximo do elemento da malha
        'Lx': 2.0,                  # Dimensão total do domínio na direção x
        'Ly': 2.0,                  # Dimensão total do domínio na direção y
        'Lz': 2.0,                  # Dimensão total do domínio na direção z
    },
    'pml': {
        'type': 'uniaxial', 
        'L': 1.0,                   # Largura da camada da PML
        'pml_order': 2,             # Ordem polinomial da PML
        'R': 1E-4                   # Coeficiente de reflexão na interface da PML
    },
    'source': {
        'type': 'gaussian',         # Tipo de fonte: 'gaussian' ou 'derivative'
        'A0': 1.0,                  # Amplitude do pulso gaussiano
        'sigma_r': 0.25,            # Desvio padrão espacial do pulso gaussiano
        'sigma_t': 1.0,             # Desvio padrão temporal do pulso gaussiano
        'x0': (1.0, 0.0, 0.0),      # Posição x do centro do pulso gaussiano
        't0': 0.0                   # Tempo de pico do pulso gaussiano
    }
}


def spatial_discretization(problem):
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
    sa = UniaxialPml(problem)
    
    # 1. Criar a malha cúbica com dimensão L com Gmsh
    gmsh.initialize()
    mesh_cubeK96_upml(problem)
    EToV = get_EToV(dim=problem['DIM'], index_based=0)
    VX, VY, VZ = extract_VX_VY_VZ(get_nodes_data(dim=problem['DIM']))
    # gmsh.fltk.run()
    gmsh.finalize()

    # 2. Criar o objeto Mesh3D 
    mesh = Mesh3D(vx=VX, vy=VY, vz=VZ, EToV=EToV, boundary_label=problem['dg']['bc'])
    # mesh.plot_mesh(title="mesh_cubeK96.msh", show_vertices=True, alpha=0.15)
    print(f"\nMalha criada com {mesh.number_of_vertices()} vértices e {mesh.number_of_elements()} elementos.")

    # 3. Definir a discretização espacial usando DG3D
    sp = Maxwell3D(
        n_order=problem['dg']['n_order'],
        mesh=mesh,
        fluxType=problem['dg']['flux_type'],
        pml_design=problem)
    
    print(f"\n🔎 Discretização espacial criada com ordem {sp.n_order}, {sp.mesh.number_of_elements()} elementos e {sp.number_of_nodes_per_element()} pontos por elemento.")

    return sa, sp


def test_with_L2_error(problem) -> None:
    """
    Testa a solução numérica comparando com a solução analítica.
    (Versão corrigida para lidar com a estrutura de dados 2D do DG-FEM)
    """
    sa, sp = spatial_discretization(problem)

    driver = MaxwellDriver(sp, CFL=problem['dg']['cfl'])
    driver['Ez'][:] = sa.gaussian_pulse_3d(sp.x, sp.y, sp.z, t=0.0)

    analytical_fields = {'Ez': sa.gaussian_pulse_3d}
    error_data = {'time': [], 'L2_error': {key: [] for key in analytical_fields}}
    
    # --- INÍCIO DA MODIFICAÇÃO CORRIGIDA DA SONDA ---    
    probe = (1.0, 0.0, 0.0)
    
    # MUDANÇA 1: A estrutura de sp.x é (nós_por_elemento, num_elementos). 
    # O cálculo da distância retorna um array com a mesma forma.
    distances = np.sqrt((sp.x - probe[0])**2 + (sp.y - probe[1])**2 + (sp.z - probe[2])**2)
    
    # MUDANÇA 2: Encontramos o índice "plano" do nó mais próximo.
    flat_probe_index = np.argmin(distances)
    
    # MUDANÇA 3: Usamos unravel_index para converter o índice plano em índices 2D 
    # (índice_do_nó_local, índice_do_elemento) que correspondem à forma de sp.x.
    node_idx, elem_idx = np.unravel_index(flat_probe_index, sp.x.shape)
    
    # MUDANÇA 4: Inicializamos a estrutura de dados da sonda
    probe_data = {'time': [], 'Ez_probe': []}
    print(f"\n🔎 Sonda pontual inserida no elemento {elem_idx}, nó local {node_idx}.")
    print(f"Coordenadas do ponto da sonda: {probe}")
    print(f"Coordenadas do nó mais próximo: ({sp.x[node_idx, elem_idx]:.3f}, {sp.y[node_idx, elem_idx]:.3f}, {sp.z[node_idx, elem_idx]:.3f})")

    # --- FIM DA MODIFICAÇÃO CORRIGIDA DA SONDA ---
    t_dg = 0.0
    N_steps = int(sa.t_final / driver.dt)
    print(f"\n🚀 Marchando no tempo com passo dt = {driver.dt:.3f}s...")
    print(f"⏳ Tempo final: {sa.t_final:.3f}s.")
    print(f"⏳ Número de passos de tempo: {N_steps}.")

    for _ in range(N_steps):
        error_data['time'].append(t_dg)
        for field_name, analytical_fn in analytical_fields.items():
            uh = driver[field_name]
            ua = analytical_fn(sp.x, sp.y, sp.z, t_dg)
            l2_error = sp.compute_L2_error(uh, ua)
            error_data['L2_error'][field_name].append(l2_error)            
        
        probe_data['time'].append(t_dg)        
        probe_value = driver['Ez'][node_idx, elem_idx]
        probe_data['Ez_probe'].append(probe_value)        
        driver.step()
        t_dg += driver.dt

    plt.figure(figsize=(10, 6))
    plt.plot(probe_data['time'], probe_data['Ez_probe'], label=f'Sonda em {probe}')
    plt.title("$Ez(x,y,z,t)$")
    plt.xlabel('Tempo (s)')
    plt.ylabel('Amplitude (V/m)')
    plt.grid(True)
    plt.legend()   


def main() -> None:
    """Função principal para execução do script."""
    clear_terminal()
    test_with_L2_error(PROBLEM)


if __name__ == '__main__':
    main()
    plt.show()
