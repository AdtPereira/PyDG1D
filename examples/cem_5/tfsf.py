#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# pylint: disable=E0611, C0413, C0103, W0401, W0614, W0622
"""
══════════════════════════════════════════════════════════════════════════

ESTE SCRIPT UTILIZA CAMINHOS ABSOLUTOS E NÃO ALTERA O DIRETÓRIO DE TRABALHO

EXECUÇÃO:
cd C:\\git\\PyDG1D
conda activate pyDG1D
python examples\\cem_5\\tfsf.py

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
from maxwell.spectralAnalyzer import PlaneWaveExcitation

from mesher.create_mesh import *
from mesher.plot_mesh import *


PROBLEM = {
    'DIM': 3,
    'description': 'Teste da formulação TF/SF do esquema DGTD tridimensional.',
    'name': 'tfsf',
    'folder': 'cem_5',
    'dg': {
        'n_order': 2,               # Ordem de interpolação polinomial
        'flux_type': 'Upwind',      # 'Upwind' or 'Centered'
        'cfl': 1.0,                 # Número de Courant-Friedrichs-Lewy
        'bc': "SMA",                # Condição de contorno: 'PEC', 'SMA'  or 'Periodic
        't_final': 8.0,             # Tempo final da simulação
    },
    'domain': {
        'type': 'cubic',            # Tipo de domínio: 'rectangle' ou 'cubic'
        'h': 4.0,                   # Tamanho máximo do elemento da malha
        'xa': 1.0,                  # Semi-lados do retângulo interno (TFZ)
        'x0': 2.0,                  # Semi-lados do retângulo intermediário (SFZ)
        'Lx': 3.0,                  # Semi-lados do retângulo externo (domínio total)
        'Ly': 3.0,                  # Dimensão total do domínio na direção y
        'Lz': 3.0,                  # Dimensão total do domínio na direção z
        'GID_TFZ': 1,               # Grupo físico para a Total Field Zone (TFZ)
    },
    'pml': {
        'type': 'uniaxial', 
        'L': 1.0,                   # Largura da camada da PML
        'pml_order': 2,             # Ordem polinomial da PML
        'R': 1E-4                   # Coeficiente de reflexão na interface da PML
    },
    'source': {
        'type': 'plane-wave',       # Tipo de fonte: 'gaussian' ou 'plane-wave'
        'E0': 1.0,                  # Amplitude do campo elétrico da fonte
        'wavelength': 1.0,          # Comprimento de onda da fonte
        'c': 1.0,                   # Velocidade da luz no meio (Normalizada)
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
    # Definição da classe PlaneWaveExcitation
    sa = PlaneWaveExcitation(problem)    

    # 1. Criar a malha cúbica com dimensão L com Gmsh
    gmsh.initialize()
    mesh_cubeK252_tfsf(problem)
    EToV = get_EToV(dim=problem['DIM'], index_based=0)
    VX, VY, VZ = extract_VX_VY_VZ(get_nodes_data(dim=problem['DIM']))
    # gmsh.fltk.run()
    gmsh.finalize()

    # 2. Criar o objeto Mesh3D 
    mesh = Mesh3D(vx=VX, vy=VY, vz=VZ, EToV=EToV, boundary_label=problem['dg']['bc'])
    print(f"\nMalha criada com {mesh.number_of_vertices()} vértices e {mesh.number_of_elements()} elementos.")    

    # 3. Plotar a malha com os ids dos pontos de colocação
    print("\n🔎 Plotando a malha cúbica K252...")
    # mesh.plot_mesh(title="mesh_cubeK254.msh", show_vertices=True, alpha=0.15)

    # 4. Criar grupos físicos na malha
    print("\n🔎 Criando grupos físicos na malha cúbica K252...")
    # Criar grupo físico cúbico TFZ de dimensão 2*xa
    # GID_TFZ = 1 -> Grupo 'TFZ'
    mesh.box_physical_group(
        group_tag=problem['domain']['GID_TFZ'],
        group_name="TFZ",
        half_dim=(sa.xa, sa.xa, sa.xa),
        center=(0.0, 0.0, 0.0),
        group_info=True
    )

    # 5. Definir a discretização espacial usando DG3D
    sp = Maxwell3D(
        n_order=problem['dg']['n_order'],
        mesh=mesh,
        fluxType=problem['dg']['flux_type'])

    # 6. Criar os perfis de condutividade da PML
    sp.sgm_x, sp.sgm_y, sp.sgm_z = sp.upml_sigma_fields(x=sp.x, y=sp.y, z=sp.z, problem=problem)

    # 7. Configura os IDs dos grupos para que 'fieldsOnInterface' saiba qual mapa usar
    # Criar a região de interface para formulação TF/SF entre os grupos 0 e 1
    sp.tfz_group_id = problem['domain']['GID_TFZ']
    sp.sf_group_id = problem['domain']['GID_PML']
    sp.buildInterfaceMaps()

    print("\n🔎 Verificando graficamente a interface TF/SF...")
    # sa.plot_interface_nodes(sp, title="Nós da Interface TF/SF (Lado do Campo Espalhado)")

    # --- INJEÇÃO DA FUNÇÃO DA FONTE ---
    sp.set_incident_wave_function(sa.get_incident_fields)   
    
    print(f"\n🔎 Discretização espacial criada com ordem {sp.n_order}, {sp.mesh.number_of_elements()} elementos e {sp.number_of_nodes_per_element()} pontos por elemento.")

    return sa, sp


def plot_probe_comparison(probe_data, probe_coords):
    """
    Plota a comparação temporal entre a solução numérica e analítica em um ponto.

    Parâmetros
    ----------
    probe_data : dict
        Dicionário com as séries temporais dos campos ('time', 'Ezh', 'Ez').
    probe_coords : tuple
        As coordenadas (x,y,z) do ponto da sonda, para o título do gráfico.
    """
    fig, ax = plt.subplots(figsize=(12, 7), dpi=100)

    # Plot da solução analítica (linha tracejada preta)
    ax.plot(probe_data['time'], probe_data['Ez'], 'k--', label='Solução Analítica ($E_z$)', linewidth=2)
    
    # Plot da solução numérica (linha azul com marcadores)
    ax.plot(probe_data['time'], probe_data['Ezh'], 'b-o', label='Solução Numérica ($E_{zh}$)', markersize=3, markeredgecolor='k', markerfacecolor='b')

    # Configurações do gráfico
    title_coords = f"({probe_coords[0]:.2f}, {probe_coords[1]:.2f}, {probe_coords[2]:.2f})"
    ax.set_title(fr'$E_z(x,y,z,t)$ @ {title_coords}', fontsize=12)
    ax.set_xlabel('Tempo (s)', fontsize=10)
    ax.set_ylabel('Amplitude (V/m)', fontsize=10)
    ax.legend(fontsize=10)
    ax.grid(True, which='both', linestyle='--', linewidth=0.5)
    ax.margins(x=0.01, y=0.05)
    plt.tight_layout()


def test_point_probe(problem, sa, sp) -> None:
    """
    Testa a solução numérica comparando com a solução analítica.
    (Versão corrigida para lidar com a estrutura de dados 2D do DG-FEM)
    """
    driver = MaxwellDriver(sp, CFL=problem['dg']['cfl']) 
    # driver['Ez'][:] = sa.incident_Ey(sp.x, sp.y, sp.z, t=0.0)
    # driver['Hz'][:] = sa.incident_Hz(sp.x, sp.y, sp.z, t=0.0)

    analytical_fields = {'Ez': sa.incident_Ey, 'Hz': sa.incident_Hz}
    error_data = {'time': [], 'L2_error': {key: [] for key in analytical_fields}}
    
    # --- INÍCIO DA SONDA ---    
    probe = (0.2, 0.0, 0.2)
    distances = np.sqrt((sp.x - probe[0])**2 + (sp.y - probe[1])**2 + (sp.z - probe[2])**2)
    flat_probe_index = np.argmin(distances)
    node_idx, elem_idx = np.unravel_index(flat_probe_index, sp.x.shape)

    # Inicializar uma estrutura de dados da sonda
    probe_data = {'time': [], 'Ezh': [], 'Hzh': [], 'Ez': [], 'Hz': []}
    print(f"\n🔎 Sonda pontual inserida no elemento {elem_idx}, nó local {node_idx}.")
    print(f"Coordenadas do ponto da sonda: {probe}")
    print(f"Coordenadas do nó mais próximo: ({sp.x[node_idx, elem_idx]:.3f}, {sp.y[node_idx, elem_idx]:.3f}, {sp.z[node_idx, elem_idx]:.3f})")

    # Coordenadas numéricas do ponto da sonda
    probe = (sp.x[node_idx, elem_idx], sp.y[node_idx, elem_idx], sp.z[node_idx, elem_idx])

    # --- INICIO DO SOLVER DG-FEM ---
    t_dg = 0.0
    N_steps = int(sa.t_final / driver.dt)
    print(f"\n🚀 Marchando no tempo com passo dt = {driver.dt:.3f}s...")
    print(f"⏳ Tempo final: {sa.t_final:.3f}s.")
    print(f"⏳ Número de passos de tempo: {N_steps}.")

    # Remova a chamada redundante para fieldsOnInterface daqui
    # corrections = sp.fieldsOnInterface(t_dg) 
    # sa.plot_incident_fields(sp, t_dg)

    for _ in range(N_steps):
        # 1. Adiciona o tempo à lista de dados da sonda no início da iteração.
        #    Isso garante que para cada conjunto de dados, haverá um ponto no tempo correspondente.
        probe_data['time'].append(t_dg)
        
        # 2. Coleta os dados do erro L2 para este instante de tempo
        for field_name, analytical_fn in analytical_fields.items():
            uh = driver[field_name]
            ua = analytical_fn(sp.x, sp.y, sp.z, t_dg)
            l2_error = sp.compute_L2_error(uh, ua)
            error_data['L2_error'][field_name].append(l2_error)            
        
        # 3. Coleta os dados numéricos e analíticos da sonda
        probe_data['Ezh'].append(driver['Ez'][node_idx, elem_idx])        
        probe_data['Hzh'].append(driver['Hz'][node_idx, elem_idx])        
        probe_data['Ez'].append(sa.incident_Ey(probe[0], probe[1], probe[2], t_dg))
        probe_data['Hz'].append(sa.incident_Hz(probe[0], probe[1], probe[2], t_dg))

        # 4. Avança a simulação no tempo
        driver.step_at_time(t_dg)
        t_dg += driver.dt

    # --- CHAMADA PARA O NOVO GRÁFICO COMPARATIVO ---
    plot_probe_comparison(probe_data, probe)


def verify_data(sa, sp):
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

    # Verificação dos dados da malha
    assert np.array_equal(sp.mesh.EToV, np.load("C:\\git\\PyDG1D\\examplesData\\outputs\\cem_5\\fpb\\EToV.npy")), "Erro: EToV."
    assert np.array_equal(sp.mesh.vx, np.load("C:\\git\\PyDG1D\\examplesData\\outputs\\cem_5\\fpb\\VX.npy")), "Erro: VX."
    assert np.array_equal(sp.mesh.vy, np.load("C:\\git\\PyDG1D\\examplesData\\outputs\\cem_5\\fpb\\VY.npy")), "Erro: VY."
    assert np.array_equal(sp.mesh.vz, np.load("C:\\git\\PyDG1D\\examplesData\\outputs\\cem_5\\fpb\\VZ.npy")), "Erro: VZ."
    assert np.array_equal(sp.EToE, np.load("C:\\git\\PyDG1D\\examplesData\\outputs\\cem_5\\fpb\\EToE.npy")), "Erro: EToE."
    assert np.array_equal(sp.EToF, np.load("C:\\git\\PyDG1D\\examplesData\\outputs\\cem_5\\fpb\\EToF.npy")), "Erro: EToF."
    print("\n🚀 Verificação dos dados da malha concluída com sucesso!")

    # Verificação dos dados da PML
    assert np.array_equal(sp.sgm_x, np.load("C:\\git\\PyDG1D\\examplesData\\outputs\\cem_5\\fpb\\sigma_x.npy")), "Erro: sgm_x."
    assert np.array_equal(sp.sgm_y, np.load("C:\\git\\PyDG1D\\examplesData\\outputs\\cem_5\\fpb\\sigma_y.npy")), "Erro: sgm_y."
    assert np.array_equal(sp.sgm_z, np.load("C:\\git\\PyDG1D\\examplesData\\outputs\\cem_5\\fpb\\sigma_z.npy")), "Erro: sgm_z."
    print("\n🚀 Verificação dos dados da PML concluída com sucesso!")

    # Verificação dos mapas de conexão
    assert np.array_equal(sp.vmapM, np.load("C:\\git\\PyDG1D\\examplesData\\outputs\\cem_5\\fpb\\vmapM.npy")), "Erro: vmapM."
    assert np.array_equal(sp.vmapP, np.load("C:\\git\\PyDG1D\\examplesData\\outputs\\cem_5\\fpb\\vmapP.npy")), "Erro: vmapP."
    assert np.array_equal(sp.vmapB, np.load("C:\\git\\PyDG1D\\examplesData\\outputs\\cem_5\\fpb\\vmapB.npy")), "Erro: vmapB."
    assert np.array_equal(sp.mapB, np.load("C:\\git\\PyDG1D\\examplesData\\outputs\\cem_5\\fpb\\mapB.npy")), "Erro: mapB."
    print("\n🚀 Verificação dos mapas de conexão concluída com sucesso!")

    print(f"\nMáscara de região total, TFZ: Dim: {sp.vmapI_G1.shape} \n{sp.vmapI_G1}")

def main() -> None:
    """Função principal para execução do script."""
    clear_terminal()
    sa, sp = spatial_discretization(PROBLEM)  
    verify_data(sa, sp)

    # test_point_probe(PROBLEM, sa, sp)
    # test_driver_time_march(PROBLEM)


if __name__ == '__main__':
    main()
    # plt.show()
