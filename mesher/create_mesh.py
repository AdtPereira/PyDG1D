import gmsh
import numpy as np
from pathlib import Path
from mesher.read_mesh import *


def mesh_rectangular_domain(PROBLEM, BOUNDARY, MATERIAL, h, view_mesh=False, mesh_info=False, auto_save=True):
    mesh_data = {}
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

    # Obter estrutura de dados DG-FEM
    VX, VY, VZ = extract_VX_VY_VZ(get_nodes_data(dim=PROBLEM['DIM'], BOUNDARY=BOUNDARY))
    EToV = get_EToV(dim=PROBLEM['DIM'], index_based=0)  # Convertendo para zero-based index

    # Create mesh Structure Data from gmsh
    mesh_data['VX'] = VX
    mesh_data['VY'] = VY
    mesh_data['EToV'] = EToV

    # Visualização da malha
    if view_mesh:
        gmsh.fltk.run()
    
    # Criação da pasta de saída
    if auto_save:
        INPUTS = (Path(__file__).parent.parent / 'examplesData' / 'inputs' / PROBLEM['folder_name']).resolve()
        INPUTS.mkdir(parents=True, exist_ok=True)
        file_path = INPUTS / f"{PROBLEM['name']}.msh"
        print(f"\nMalha salva em {file_path}")
        gmsh.write(str(file_path))
        # basic_info()

    # Exibir informações da malha
    if mesh_info:
        print("\n🌐 Informações da malha:")
        print(f"🌐 VX:\n {mesh_data['VX']}")
        print(f"🌐 VY:\n {mesh_data['VY']}")
        print(f"🌐 VZ:\n {mesh_data['VZ']}")
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


def mesh_single_conductor_domain(problem, h, view_mesh=False, auto_save=True):
    mesh_data = {}
    TRIANGLE_TYPE = 2

    # Inicializar o Gmsh
    gmsh.initialize()
    gmsh.model.add("single_conductor_domain")
    factory = gmsh.model.occ

    # ----------------------------
    # Parâmetros da geometria
    # ----------------------------
    rc = problem['domain']['rc']    # raio do disco central
    x0 = problem['pml']['x0']       # semi-lados do retângulo interno
    y0 = problem['pml']['y0']       # semi-lados do retângulo interno
    L = problem['pml']['L']         # camada externa da PML

    # ----------------------------
    # Geometria base
    # ----------------------------
    outer = factory.addRectangle(-(x0 + L), -(y0 + L), 0, 2*(x0 + L), 2*(y0 + L))
    inner = factory.addRectangle(-x0, -y0, 0, 2*x0, 2*y0)
    disk  = factory.addDisk(0, 0, 0, rc, rc)

    # ----------------------------
    # Operações booleanas (sem sobreposição)
    # ----------------------------
    region_inner, _ = factory.cut([(2, inner)], [(2, disk)], removeObject=True, removeTool=True)
    region_outer, _ = factory.cut([(2, outer)], [(2, inner)], removeObject=True, removeTool=False)
    region_core  = [(2, disk)]

    # ----------------------------
    # Malha
    # ----------------------------
    factory.synchronize()
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
        INPUTS = (Path(__file__).parent.parent / 'examplesData' / 'inputs' / problem['folder_name']).resolve()
        INPUTS.mkdir(parents=True, exist_ok=True)
        file_path = INPUTS / f"{problem['name']}.msh"
        print(f"\nMalha salva em {file_path}")
        gmsh.write(str(file_path))
        # basic_info()

    # Create mesh Structure Data from gmsh
    mesh_data['VX'] = VX
    mesh_data['VY'] = VY
    mesh_data['EToV'] = EToV

    gmsh.finalize()
    print(f"\n🌐 Malha criada com {len(VX)} nós e {len(EToV)} elementos triangulares.\n")
    return mesh_data


def mesh_cubeK6():
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


def mesh_cubeK24(problem):
    # --- Parâmetros Geométricos e de Malha ---
    Lx = problem['domain']['Lx']    # Dimensão total do domínio na direção x
    Ly = problem['domain']['Ly']    # Dimensão total do domínio na direção y
    Lz = problem['domain']['Lz']    # Dimensão total do domínio na direção z

    gmsh.model.add("cubeK24")

    # 1. Cria o cubo da mesma forma que antes
    box = gmsh.model.occ.addBox(0, 0, 0, Lx, Ly, Lz)
    gmsh.model.occ.synchronize()

    # 5. Gera a malha 3D
    gmsh.option.setNumber("Mesh.MeshSizeMax", value=1)
    gmsh.option.setNumber("Mesh.MeshSizeMin", value=1)
    gmsh.model.mesh.generate(3)


def mesh_cubeK96_upml(problem):
    # --- Parâmetros Geométricos e de Malha ---
    L = problem['pml']['L']         # Largura da camada da PML
    h = problem['domain']['h']      # Tamanho máximo do elemento da malha
    Lx = problem['domain']['Lx']    # Dimensão total do domínio na direção x
    x0 = Lx - L                     # semi-lados do retângulo interno - Domínio Físico

    gmsh.model.add("cubeK96_upml")
    factory = gmsh.model.occ

    # --- Criação da Geometria ---
    tfz = factory.addBox(-x0, -x0, -x0, 2*x0, 2*x0, 2*x0)
    pml = factory.addBox(-Lx, -Lx, -Lx, 2*Lx, 2*Lx, 2*Lx)

    pml_domain, _ = factory.cut([(3, pml)], [(3, tfz)], removeTool=False)

    factory.synchronize()

    # --- Definição dos Grupos Físicos ---
    # Extrai a tag numérica de cada volume criado.
    # O resultado das operações de corte pode ter mais de um volume, mas aqui sabemos que é apenas um.
    pmlTag = pml_domain[0][1]
    # A tag do domínio total-field é a original, que geralmente é 1
    tfzTag = tfz

    # Atribui cada volume a um grupo físico com um nome e uma tag numérica única.
    # A dimensão é 3 para volumes.
    gmsh.model.addPhysicalGroup(3, [tfzTag], tag=301, name="TFZ")
    gmsh.model.addPhysicalGroup(3, [pmlTag], tag=303, name="PML")

    # --- Malha ---
    gmsh.option.setNumber("Mesh.MeshSizeMax", h)
    gmsh.option.setNumber("Mesh.MeshSizeMin", h)
    gmsh.model.mesh.generate(3)


def mesh_cubeK252_tfsf(problem):
    # --- Parâmetros Geométricos e de Malha ---
    h = problem['domain']['h']      # Tamanho máximo do elemento da malha
    xa = problem['domain']['xa']    # Semi-lado do retângulo interno (TFZ)
    x0 = problem['domain']['x0']    # Semi-lado do retângulo intermediário (SFZ)
    Lx = problem['domain']['Lx']    # Semi-lado do retângulo externo (domínio total)

    gmsh.model.add("cubeK168_tfsf")
    factory = gmsh.model.occ

    # --- Criação da Geometria ---
    tfz = factory.addBox(-xa, -xa, -xa, 2*xa, 2*xa, 2*xa)
    sfz = factory.addBox(-x0, -x0, -x0, 2*x0, 2*x0, 2*x0)
    pml = factory.addBox(-Lx, -Lx, -Lx, 2*Lx, 2*Lx, 2*Lx)

    # Realiza o corte do volume PML pelo volume SFZ
    pml_domain, _ = factory.cut([(3, pml)], [(3, sfz)], removeTool=False)
    sfz_domain, _ = factory.cut([(3, sfz)], [(3, tfz)], removeTool=False)

    factory.synchronize()

    # --- Definição dos Grupos Físicos ---
    # Extrai a tag numérica de cada volume criado.
    # O resultado das operações de corte pode ter mais de um volume, mas aqui sabemos que é apenas um.
    pmlTag = pml_domain[0][1]
    sfzTag = sfz_domain[0][1]
    # A tag do domínio total-field é a original, que geralmente é 1
    tfzTag = tfz

    # Atribui cada volume a um grupo físico com um nome e uma tag numérica única.
    # A dimensão é 3 para volumes.
    gmsh.model.addPhysicalGroup(3, [tfzTag], tag=301, name="TFZ")
    gmsh.model.addPhysicalGroup(3, [sfzTag], tag=302, name="SFZ")
    gmsh.model.addPhysicalGroup(3, [pmlTag], tag=303, name="PML")

    # --- Malha ---
    gmsh.option.setNumber("Mesh.MeshSizeMax", h)
    gmsh.option.setNumber("Mesh.MeshSizeMin", h)
    gmsh.model.mesh.generate(3)