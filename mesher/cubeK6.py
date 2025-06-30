import os
import gmsh

from .read_mesh import *
from .create_mesh import mesh_cubeK6    

if __name__ == "__main__":
    # Clear terminal
    os.system('cls' if os.name == 'nt' else 'clear')
    DIM = 3

    # 1. Inicializa o GMSH
    gmsh.initialize()

    # 2. Define o modelo
    mesh_cubeK6()

    # 3, Imprime informações da malha
    basic_info(dim=DIM)
    # complete_info()

    # 4. Obtém a conectividade dos elementos
    EToV = get_EToV(dim=DIM, index_based=0)  # Convertendo para zero-based index
    print(f"\nMatriz EToV de Conectividade (dim={DIM}): \n{EToV}")

    # 5. Obtém os dados dos nós
    nodes_data = get_nodes_data(dim=DIM)
    print(f"\nDados dos nós (dim={DIM}):")
    for node, data in nodes_data.items():
        print(f"Node {node}: xg={data['xg']}")

    # 6. Obtém as estruturas de dados do DG-FEM
    VX, VY, VZ = extract_VX_VY_VZ(nodes_data)
    print(f"\nCoordenadas dos nós (dim={DIM}):")
    print(f"VX: {VX}")
    print(f"VY: {VY}")  
    print(f"VZ: {VZ}")

    # 6. Exibe a malha e Finaliza o GMSH
    gmsh.fltk.run()
    gmsh.finalize()
