import gmsh
import os

# Limpar o terminal
os.system('cls' if os.name == 'nt' else 'clear')

# Define os parâmetros de entrada
FINITE_ELEMENT = ("Tetrahedron", 1)
BOUNDARY = [{'tag': 201, 'type': 'Dirichlet', 'value': 0.0, 'name': 'nxE=0'}]
MATERIAL = [{'tag': 301, 'name': 'free_space', 'ur': 1, 'er': 1}]

# Tamanho do elemento da malha
element, order = FINITE_ELEMENT
h = 1  

# Dimensões da cavidade retangular, metros
a, b, c = 1.0, 1.0, 1.0

# Inicializar o Gmsh
gmsh.initialize()
gmsh.model.add("rectangular_cavity")

# Criar volume retangular (caixa)
TagVolume = gmsh.model.occ.addBox(0, 0, 0, a, b, c)
print(f"\nTagVolume: {TagVolume}")
gmsh.model.occ.synchronize()

gmsh.option.setNumber("Mesh.MeshSizeMax", h)
gmsh.option.setNumber("Mesh.MeshSizeMin", h)
gmsh.model.mesh.generate(dim=3)
gmsh.model.mesh.setOrder(order)

print("\n----------- Beginning of edge mapping --------------------")
# If all you need is the list of all edges or faces in terms of their nodes:
# gmsh.model.mesh.createEdges([(3, 1)])
gmsh.model.mesh.createEdges()
edgeTags, edgeNodes = gmsh.model.mesh.getAllEdges()
print(f"Total number of edges: {len(edgeTags)}")
print(f"Total number of edgeNodes: {len(edgeNodes)}")

# Criar edge_mapping
edge_mapping = {tuple(sorted([edgeNodes[2*i], edgeNodes[2*i + 1]])): tag
                 for i, tag in enumerate(sorted(edgeTags))}
print(f"edge_mapping:\n {edge_mapping}")
print("----------- Ending of edge mapping --------------------")

print("\n----------- Beginning of Face mapping --------------------")
# If all you need is the list of all edges or faces in terms of their nodes:
# gmsh.model.mesh.createEdges([(3, 1)])
gmsh.model.mesh.createFaces()
faceTags, faceNodes = gmsh.model.mesh.getAllFaces(3)
print(f"Total number of faces: {len(faceTags)}")
print(f"Total number of faceNodes: {len(faceNodes)}")

# Criar face_mapping
face_mapping = {tuple(sorted([faceNodes[3*i], faceNodes[3*i + 1], faceNodes[3*i + 2]])): tag
                 for i, tag in enumerate(sorted(faceTags))}
print(f"face_mapping:\n {face_mapping}")
print("----------- Ending of edge mapping --------------------")

print("\n----------- Beginning of Boundary mapping --------------------")
# Obter os contornos (superfícies, dim=2) do volume
BoundaryDimTags = gmsh.model.getBoundary([(3, TagVolume)], oriented=True, recursive=False)
print(f"BoundaryDimTags: {BoundaryDimTags}")

# Exibir os TAGs das superfícies associadas a cada contorno
BoundaryTags = [Dimtags[1] for Dimtags in BoundaryDimTags]

# Definindo as superfícies de contorno de Dirichlet (dim=2)
gmsh.model.addPhysicalGroup(dim=2, tags=BoundaryTags, tag=BOUNDARY[0]['tag'], name=BOUNDARY[0]['name'])

# Adicionar grupos físicos para Dim=3 (volume)
gmsh.model.addPhysicalGroup(dim=3, tags=[TagVolume], tag=MATERIAL[0]['tag'], name=MATERIAL[0]['name'])

## ------------------##
## `get_cell_data()` ##
## ------------------##

# 1. Obter os elementos da malha
elementTypes, elementTags, nodeTags = gmsh.model.mesh.getElements(dim=3)
print(f"elemTypes with dim 3: {elementTypes}")
print(f"elemTags with dim 3: {elementTags}")
print(f"len(elemTags[0]): {len(elementTags[0])}")
print(f"elemNodeTags: {nodeTags}")
print(f"len(nodeTags[0]): {len(nodeTags[0])}")
conn_dict = {}

# 2. Obter as entidades geometricas do modelo
print(f"Geometrical Entities with dim 0 (Points): {gmsh.model.getEntities(dim=0)}")
print(f"Geometrical Entities with dim 1 (Edges): {gmsh.model.getEntities(dim=1)}")
print(f"Geometrical Entities with dim 2 (Faces): {gmsh.model.getEntities(dim=2)}")
print(f"Geometrical Entities with dim 3 (Regions): {gmsh.model.getEntities(dim=3)}")
print("----------- Ending of edge mapping --------------------")

print("\n----------- Beginning Material Physical groups with dim3 --------------------")
# 3. Criar o dicionário mesh_data['cell]
cell_data = {}

for material in MATERIAL:
    # Obter as entidades físicas (grupo físico) associadas ao material
    MaterialEntitiesTags = gmsh.model.getEntitiesForPhysicalGroup(dim=3, tag=material['tag'])
    print(f"Entities of {material['name']} with dim3: {MaterialEntitiesTags}")
    
    for EntityTag in MaterialEntitiesTags:
        elementTypes, elementTags, nodeTags = gmsh.model.mesh.getElements(dim=3, tag=EntityTag)
        print(f"mesh_elements of entity (dim, tag) = (3, {EntityTag}): elemTags = {elementTags}")
        print(f"mesh_elements of entity (dim, tag) = (3, {EntityTag}): nodeTags = {nodeTags}")
        
        for elemType, elemTag, elemNode in zip(elementTypes, elementTags, nodeTags):            
            # Obter as propriedades do elemento
            _, _, _, nodes_per_element, _, _ = gmsh.model.mesh.getElementProperties(elemType)

            # Número de elementos
            N_tet = len(elemNode) // nodes_per_element

            # Criar dicionário associando cada elemento à sua conectividade de nós
            for i, Tag in enumerate(elemTag):                
                # Obter a conectividade do elemento
                # conn_node = list(sorted(elemNodes[0][4 * i : 4 * (i + 1)])) 
                conn_node = elemNode[nodes_per_element * i: nodes_per_element * (i + 1)].tolist()
                conn_std = sorted(conn_node)

                # Tetrahedron element
                if nodes_per_element == 4:
                    conn_edge = [
                        edge_mapping[(conn_std[0], conn_std[1])],  # e1: 1 -> 2
                        edge_mapping[(conn_std[0], conn_std[2])],  # e2: 1 -> 3
                        edge_mapping[(conn_std[0], conn_std[3])],  # e3: 1 -> 4
                        edge_mapping[(conn_std[1], conn_std[2])],  # e4: 2 -> 3
                        edge_mapping[(conn_std[1], conn_std[3])],  # e5: 2 -> 4
                        edge_mapping[(conn_std[2], conn_std[3])]   # e6: 3 -> 4
                    ]

                    # Conectividade de faces (cada face tem 3 nós)
                    conn_face = [
                        sorted([conn_node[0], conn_node[1], conn_node[2]]),  # Face 1
                        sorted([conn_node[0], conn_node[1], conn_node[3]]),  # Face 2
                        sorted([conn_node[0], conn_node[2], conn_node[3]]),  # Face 3
                        sorted([conn_node[1], conn_node[2], conn_node[3]])   # Face 4
                    ]
                
                # Adicionar ao dicionário de células
                cell_data[Tag] = {
                    'tag': Tag,
                    'conn': conn_node,
                    'conn_sorted': conn_std,
                    'conn_edge': conn_edge,
                    'conn_face': conn_face,
                    # 'geo': {'centroid': None, 'dim': None},
                    # 'contour': {'type': None, 'conn_contour': None},
                    'material': material['tag']}

print(f"cell_data for dim3: {cell_data}")
print("----------- Ending Material Physical groups --------------------")

print("\n----------- Beginning Boundary Physical groups with dim2 --------------------")
# 4. Criar um mapa entre cada entidade física (grupo físico) e os elementos correspondentes
boundary_data = {}

for bc in BOUNDARY:
    # Obter as entidades físicas (grupo físico) associadas ao contorno
    BoundaryEntitiesTags = gmsh.model.getEntitiesForPhysicalGroup(dim=2, tag=bc['tag'])
    print(f"Entities of {bc['name']} with dim2: {BoundaryEntitiesTags}")
    
    for EntityTag in BoundaryEntitiesTags:
        elementTypes, elementTags, nodeTags = gmsh.model.mesh.getElements(dim=2, tag=EntityTag)
        print(f"mesh_elements of entity (dim, tag) = (2, {EntityTag}): elemTags = {elementTags}")
        print(f"mesh_elements of entity (dim, tag) = (2, {EntityTag}): nodeTags = {nodeTags}")
        
        for elemType, elemTag, elemNode in zip(elementTypes, elementTags, nodeTags):   
            # Obter as propriedades do elemento
            _, _, _, nodes_per_element, _, _ = gmsh.model.mesh.getElementProperties(elemType)

            # Número de elementos
            N_tri = len(elemNode) // nodes_per_element             

            # Criar dicionário associando cada elemento à sua conectividade de nós
            for i, Tag in enumerate(elemTag):    
                # Obter a conectividade do elemento            
                # conn_node = list(sorted(nodeTags[0][3 * i : 3 * (i + 1)]))
                conn_node = elemNode[nodes_per_element * i: nodes_per_element * (i + 1)].tolist()
                conn_std = sorted(conn_node)

                # Triangular element
                if nodes_per_element == 3:
                    conn_edge = [
                        edge_mapping[(conn_std[0], conn_std[1])],  # e1: 1 -> 2
                        edge_mapping[(conn_std[0], conn_std[2])],  # e2: 1 -> 3
                        edge_mapping[(conn_std[1], conn_std[2])],  # e3: 2 -> 3
                    ]

                # Dicionário de contorno
                boundary_data[Tag] = {
                    'tag': Tag,
                    'conn': conn_node,
                    'conn_sorted': conn_std,
                    'conn_edge': conn_edge,
                    'boundary': bc['tag']}

print(f"boundary data for dim2: {boundary_data}")
print("----------- Ending Boundary Physical groups with dim2 --------------------")

# Visualizar a malha no ambiente Gmsh (opcional)
gmsh.fltk.run()

# Finalizar o Gmsh
gmsh.finalize()
