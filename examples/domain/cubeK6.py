import os
import gmsh
from mesher.read_mesh import *

if __name__ == "__main__":
    # Clear terminal
    os.system('cls' if os.name == 'nt' else 'clear')
    DIM = 3

    # 1. Inicializa o GMSH
    gmsh.initialize()

    # 2. Define o modelo
    gmsh.model.add("cubeK6")

    # 1. Cria o cubo da mesma forma que antes
    box = gmsh.model.occ.addBox(0, 0, 0, 1, 1, 1)
    gmsh.model.occ.synchronize()

    # 2. Força uma malha estruturada (Transfinita)    
    # Obtém as "entidades de fronteira" (as 6 faces) do nosso volume.
    surfaces = gmsh.model.getBoundary([(3, box)], combined=False, oriented=False)
    
    # Extrai apenas as tags numéricas das superfícies para uso posterior
    surface_tags = [s[1] for s in surfaces]

    for s_tag in surface_tags:
        # Para cada face, define que ela será transfinita.
        gmsh.model.mesh.setTransfiniteSurface(s_tag)
        
        # Pega as arestas (curves) da face atual, também com 'oriented=False'.
        curves = gmsh.model.getBoundary([(2, s_tag)], combined=False, oriented=False)
        for c_tag in [c[1] for c in curves]:
            # Define que cada aresta terá 2 pontos (os pontos inicial e final)
            gmsh.model.mesh.setTransfiniteCurve(c_tag, 2)
            
    # Define o volume como transfinito
    gmsh.model.mesh.setTransfiniteVolume(box)
    
    # 3. Força a recombinação para criar hexaedros primeiro
    gmsh.model.mesh.setRecombine(3, box)
    
    # 4. Força a subdivisão dos hexaedros em tetraedros
    gmsh.option.setNumber("Mesh.SubdivisionAlgorithm", 1)

    # 5. Adiciona um grupo físico para as superfícies de contorno
    # A função addPhysicalGroup recebe a dimensão (2 para superfícies),
    # uma lista com as tags das entidades e uma tag para o grupo físico.
    gmsh.model.addPhysicalGroup(2, surface_tags, tag=201, name="BoundarySurfaces")

    # 6. Gera a malha 3D
    gmsh.model.mesh.generate(3)

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
