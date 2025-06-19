import gmsh
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from mesher import read_mesh


def mesh_rectangular_domain(PROBLEM, BOUNDARY, MATERIAL, h, view_mesh=False, mesh_info=False, auto_save=True):
    mesh_data = {}
    TRI_TYPE = 2

    # Dimensões do domínio retangular
    a, b = PROBLEM['Lx'], PROBLEM['Ly']

    # Inicializar o Gmsh
    gmsh.initialize()
    gmsh.model.add("rectangular_domain")

    # Criar superfície retangular
    TagSurface = gmsh.model.occ.addRectangle(-a/2, -b/2, 0, a, b)
    gmsh.model.occ.synchronize()
    gmsh.option.setNumber("Mesh.MeshSizeMin", h)
    gmsh.option.setNumber("Mesh.MeshSizeMax", h)
    gmsh.model.mesh.generate(dim=2)
    gmsh.model.mesh.setOrder(1)

    # Obter os contornos (curvas, dim=1) de cada superfície
    outDimTags = gmsh.model.getBoundary([(2, TagSurface)], oriented=True, recursive=False)

    # Exibir os TAGs das curvas associadas a cada contorno
    tagList_boundary = [Dimtags[1] for Dimtags in outDimTags]

    # Definindo as curvas de contorno de Dirichlet (dim=1)
    gmsh.model.addPhysicalGroup(dim=1, tags=tagList_boundary, tag=BOUNDARY[0]['tag'], name=BOUNDARY[0]['name'])

    # Adicionar grupos físicos para Dim=2 (superfícies)
    gmsh.model.addPhysicalGroup(dim=2, tags=[TagSurface], tag=MATERIAL[0]['tag'], name=MATERIAL[0]['name'])

    # Obter dados de nós (vértices)
    _ , node_coords, _ = gmsh.model.mesh.getNodes()
    coords = node_coords.reshape(-1, 3)

    VX = coords[:, 0]  # x-coordinates
    VY = coords[:, 1]  # y-coordinates

    # Obter elementos de dimensão 2 (triângulos)
    elem_types, _, elem_node_tags = gmsh.model.mesh.getElements(dim=2)
    
    # Considerar apenas elementos triangulares (type 2)
    if TRI_TYPE in elem_types:
        idx = np.where(elem_types == TRI_TYPE)[0][0]
        elem_nodes = elem_node_tags[idx]
        EToV = np.array(elem_nodes, dtype=int).reshape(-1, 3) - 1  # Convertendo para zero-based index
    else:
        raise ValueError("Malha não contém elementos triangulares.")
    
    if view_mesh:
        gmsh.fltk.run()
    
    # Criação da pasta de saída, se habilitado
    if auto_save:
        INPUTS = (Path(__file__).parent.parent / 'examplesData' / 'inputs' / PROBLEM['folder_name']).resolve()
        INPUTS.mkdir(parents=True, exist_ok=True)
        file_path = INPUTS / f"{PROBLEM['name']}.msh"
        print(f"\nMalha salva em {file_path}")
        gmsh.write(str(file_path))
        # basic_info()

    # Create mesh Structure Data from gmsh
    mesh_data['cell'] = read_mesh.get_cell_data(MATERIAL)
    mesh_data['nodes'] = read_mesh.get_nodes_data(BOUNDARY, problem_dim=2)
    mesh_data['VX'] = VX
    mesh_data['VY'] = VY
    mesh_data['EToV'] = EToV

    # Exibir informações da malha
    if mesh_info:
        print(f"🌐 VX: {mesh_data['VX']}")
        print(f"🌐 VY: {mesh_data['VY']}")
        print(f"🌐 EToV:\n {mesh_data['EToV']}")

    gmsh.finalize()
    print(f"\n🌐 Malha criada com {len(VX)} nós e {len(EToV)} elementos triangulares.\n")
    return mesh_data


def mesh_rectangular_pml_domain(PROBLEM, PML_DESIGN, BOUNDARY, MATERIAL, h, view_mesh=False, mesh_info=False, auto_save=True):
    mesh_data = {}
    TRIANGLE_TYPE = 2

    # Inicializar o Gmsh
    gmsh.initialize()
    gmsh.model.add("rectangular_pml")
    factory = gmsh.model.occ

    # Dimensões do domínio 
    L = PML_DESIGN['L']         # Espessura da PML
    x0 = PML_DESIGN['x0']       # Fronteira interna da PML
    xa, xb = x0, x0 + L         # Fronteiras externas da PML

    # Criar regiões absorvedoras, omega_PML
    region_a = factory.addRectangle(-xb, -xb, 0, L, L)
    region_i = factory.addRectangle(-xa, -xb, 0, 2*x0, L)
    region_b = factory.addRectangle(xa, -xb, 0, L, L)
    region_ii = factory.addRectangle(xa, -xa, 0, L, 2*x0)
    region_c = factory.addRectangle(xa, xa, 0, L, L)
    region_iii = factory.addRectangle(-xa, xa, 0, 2*x0, L)
    region_d = factory.addRectangle(-xb, xa, 0, L, L)
    region_iv = factory.addRectangle(-xb, -xa, 0, L, 2*x0)

    # Criar região do espaço livre, omega_fs
    region_fs = factory.addRectangle(-xa, -xa, 0, 2*x0, 2*x0)

    # Fragmentar todas as regiões para garantir interfaces conformais
    objectDimTags = [
        (2, region_fs),
        (2, region_i), (2, region_ii), (2, region_iii), (2, region_iv), 
        (2, region_a), (2, region_b), (2, region_c), (2, region_d)
    ]

    factory.fragment(objectDimTags, objectDimTags)
    factory.synchronize()

    # Adicionar grupos físicos para superfícies (Dim=2)	    
    TagList_surfaces = [region_fs, region_a, region_b, region_c, region_d, region_i, region_ii, region_iii, region_iv]
    for i, SurfaceList in enumerate(TagList_surfaces):
        gmsh.model.addPhysicalGroup(2, [SurfaceList], tag=MATERIAL[i]['tag'], name=MATERIAL[i]['name'])

    # Definir ordem dos elementos
    gmsh.option.setNumber("Mesh.MeshSizeMax", h)
    gmsh.option.setNumber("Mesh.MeshSizeMin", h)
    gmsh.model.mesh.generate(2)
    gmsh.model.mesh.setOrder(1)

    # Obter dados de nós (vértices)
    _ , node_coords, _ = gmsh.model.mesh.getNodes()
    coords = node_coords.reshape(-1, 3)

    VX = coords[:, 0]  # x-coordinates
    VY = coords[:, 1]  # y-coordinates

    # Obter elementos de dimensão 2 (triângulos)
    elem_types, _, elem_node_tags = gmsh.model.mesh.getElements(dim=2)
    
    # Considerar apenas elementos triangulares (type 2)
    if TRIANGLE_TYPE in elem_types:
        idx = np.where(elem_types == TRIANGLE_TYPE)[0][0]
        elem_nodes = elem_node_tags[idx]
        EToV = np.array(elem_nodes, dtype=int).reshape(-1, 3) - 1  # Convertendo para zero-based index
    else:
        raise ValueError("Malha não contém elementos triangulares.")
    
    if view_mesh:
        gmsh.fltk.run()
    
    # Criação da pasta de saída, se habilitado
    if auto_save:
        INPUTS = (Path(__file__).parent.parent / 'examplesData' / 'inputs' / PROBLEM['folder_name']).resolve()
        INPUTS.mkdir(parents=True, exist_ok=True)
        file_path = INPUTS / f"{PROBLEM['name']}.msh"
        print(f"\nMalha salva em {file_path}")
        gmsh.write(str(file_path))
        # basic_info()

    # Create mesh Structure Data from gmsh
    mesh_data['cell'] = read_mesh.get_cell_data(MATERIAL)
    mesh_data['nodes'] = read_mesh.get_nodes_data(BOUNDARY, problem_dim=2)
    mesh_data['VX'] = VX
    mesh_data['VY'] = VY
    mesh_data['EToV'] = EToV

    # Exibir informações da malha
    if mesh_info:
        print(f"🌐 VX: {mesh_data['VX']}")
        print(f"🌐 VY: {mesh_data['VY']}")
        print(f"🌐 EToV:\n {mesh_data['EToV']}")

    gmsh.finalize()
    print(f"\n🌐 Malha criada com {len(VX)} nós e {len(EToV)} elementos triangulares.\n")
    return mesh_data


def structured_data(mesh_data):
    # Extraindo as coordenadas globais dos nós (x, y)
    try:
        xg = [node['xg'][0] for node in mesh_data['nodes'].values()]
        yg = [node['xg'][1] for node in mesh_data['nodes'].values()]
    except KeyError as e:
        raise ValueError(f"Erro ao acessar as coordenadas globais: {e}")

    # Extraindo a matriz de conectividade
    try:
        conn = [cell['conn'] for cell in mesh_data['cell'].values()]
        # Ajusta índice para 0-based
        conn_py = [[node - 1 for node in nodes[:3]] for nodes in conn]  
    except KeyError as e:
        raise ValueError(f"Erro ao acessar a conectividade: {e}")
    
    return xg, yg, conn_py


def plot_triangular_mesh(INFO_GRAPH, mesh_data, nodes_index_based=0):
    # Dados do gráfico
    show_cell = INFO_GRAPH['cell']
    show_nodes = INFO_GRAPH['nodes']
    show_edges = INFO_GRAPH['edges']
    show_edges_numb = INFO_GRAPH['edges_numb']
    filepath = INFO_GRAPH['filepath']

    # Estruturando os dados da malha
    nodes_data = mesh_data['nodes']
    
    # Extraindo as coordenadas globais dos nós (x, y) e a matriz de conectividade
    xg, yg, conn_py = structured_data(mesh_data)

    # Plotando a malha de elementos finitos
    plt.figure(figsize=(8, 6))
    plt.triplot(xg, yg, conn_py, color='gray')  
    
    # Plotando as arestas dos elementos
    if show_edges or show_edges_numb:
        for key, edge in mesh_data['edges'].items():
            # Coordenadas dos nós inicial e final
            x0 = nodes_data[edge['conn'][0]]['xg'][0]
            y0 = nodes_data[edge['conn'][0]]['xg'][1]
            x1 = nodes_data[edge['conn'][1]]['xg'][0]
            y1 = nodes_data[edge['conn'][1]]['xg'][1]
            
            # Ponto médio da aresta
            x_mid, y_mid = (x0 + x1) / 2, (y0 + y1) / 2
            
            # Vetor da seta (a partir do ponto médio)
            dx, dy = (x1 - x0) * 0.2, (y1 - y0) * 0.2  
                
            # Adicionando uma seta no meio da aresta
            if show_edges:
                plt.arrow(x_mid, y_mid, dx, dy, head_width=0.015, head_length=0.05,
                                fc='blue', ec='blue', length_includes_head=True)

            # Adicionando os números das arestas
            if show_edges_numb:
                plt.scatter(x_mid, y_mid, marker='s', color='white', 
                                edgecolor='black', s=120, zorder=1)                
                plt.text(x_mid, y_mid, key, color='blue', fontsize=6, 
                                ha='center', va='center')

    # Adicionando nós
    if show_nodes: 
        for key, node in nodes_data.items():
            x, y = node['xg'][0], node['xg'][1]
            plt.scatter(x, y, color='white', edgecolor='black', s=180)
            if nodes_index_based == 1:
                plt.text(x, y, str(key), color='red', fontsize=6, ha='center', va='center')
            elif nodes_index_based == 0:
                plt.text(x, y, str(key - 1), color='red', fontsize=6, ha='center', va='center')
            else:
                raise ValueError("nodes_index_based deve ser 0 ou 1.")  
    else:
        plt.scatter(xg, yg, color='black', s=1, zorder=3)

    # Adicionando elementos
    if show_cell:
        for key, cell in mesh_data['cell'].items():
            xc = np.mean([nodes_data[node]['xg'][0] for node in cell['conn']])
            yc = np.mean([nodes_data[node]['xg'][1] for node in cell['conn']])
            plt.text(xc, yc, str(key), fontweight='bold',
                        color='black', fontsize=6, ha='center', va='center')
                
    # Ajustando rótulos e layout
    plt.xlabel(r'$x$')
    plt.ylabel(r'$y$')
    plt.axis('equal')
    plt.tight_layout()
    
    # Salvando o arquivo no formato SVG
    plt.savefig(filepath, format="svg")
    plt.show()
    print(f"Arquivo gráfico da malha salvo em: {filepath}")