import gmsh
import os

# Limpar o terminal
os.system('cls' if os.name == 'nt' else 'clear')

BOUNDARY = [{'tag': 101, 'type': 'PML', 'value': None, 'name': 'inner_pml_domain'},
            {'tag': 102, 'type': 'PML', 'value': None, 'name': 'outer_pml_domain'}
]

MATERIAL = [{'tag': 201, 'name': 'free_space', 'relative_magnetic_permeability': 1, 'relative_electric_permittivity': 1},
            {'tag': 202, 'name': 'pml', 'relative_magnetic_permeability': 1, 'relative_electric_permittivity': 1}
]

WAVELENGTH = 1.0  # Comprimento de onda em metros
GEOMETRY = {'a': {'h': WAVELENGTH, 'L': WAVELENGTH, 'x0': WAVELENGTH}}

# Inicializar o Gmsh
gmsh.initialize()
gmsh.model.add("rectangular_pml")
factory = gmsh.model.occ

# Dimensões do domínio 
h = GEOMETRY['a']['h']    # Tamanho do elemento
L = GEOMETRY['a']['L']    # Lado do retângulo externo
x0 = GEOMETRY['a']['x0']  # Lado do retângulo interno
y0 = x0

# Criar as superfícies retangulares externa e interna
rect_outer = factory.addRectangle(-(x0 + L), -(y0 + L), 0, 2 * (x0 + L), 2 * (y0 + L))
rect_inner = factory.addRectangle(-x0, -y0, 0, 2 * x0, 2 * y0)

# Free Space
outDimTags_fs = [(2, rect_inner)]
print(f"Tag da superfície 'free_space': {outDimTags_fs}")

# Subtrair o retângulo interno do retângulo externo para criar a PML
# O retângulo interno (ferramenta) não é removido (removeTool=False)
outDimTags_pml, _ = factory.cut([(2, rect_outer)], outDimTags_fs, removeTool=False)
print(f"Tag da superfície 'pml': {outDimTags_pml}")

# Sincronizar após a operação de corte
factory.synchronize()

# --- INÍCIO DA CORREÇÃO ---

# Obter as curvas de fronteira de cada superfície
boundary_pml_all = gmsh.model.getBoundary(outDimTags_pml, oriented=False, recursive=False)
boundary_inner = gmsh.model.getBoundary(outDimTags_fs, oriented=False, recursive=False)

print(f"Curvas da fronteira interna (rect_inner): {boundary_inner}")
print(f"Todas as curvas da fronteira da PML: {boundary_pml_all}")

# Extrair apenas os tags numéricos para facilitar a manipulação
tagList_inner = [tag[1] for tag in boundary_inner]
tagList_pml_all = [tag[1] for tag in boundary_pml_all]

# Para obter as curvas da fronteira externa, fazemos a diferença entre
# o conjunto de todas as curvas da PML e o conjunto das curvas internas.
tagList_outer = list(set(tagList_pml_all) - set(tagList_inner))

print(f"Tags das curvas da fronteira interna (isoladas): {tagList_inner}")
print(f"Tags das curvas da fronteira externa (isoladas): {tagList_outer}")

# Adicionar grupos físicos para curvas (Dim=1)
print("\nGrupos físicos de Curvas (Dim=1):")
# Fronteira Interna
gmsh.model.addPhysicalGroup(1, tagList_inner, tag=BOUNDARY[0]['tag'], name=BOUNDARY[0]['name'])
print(f"  Tag {BOUNDARY[0]['tag']} ('{BOUNDARY[0]['name']}'): {tagList_inner}")
# Fronteira Externa
gmsh.model.addPhysicalGroup(1, tagList_outer, tag=BOUNDARY[1]['tag'], name=BOUNDARY[1]['name'])
print(f"  Tag {BOUNDARY[1]['tag']} ('{BOUNDARY[1]['name']}'): {tagList_outer}")

# Adicionar grupos físicos para superfícies (Dim=2)
print("\nGrupos físicos de Superfícies (Dim=2):")
# Superfície 'free_space'
tag_fs_surface = [dimTag[1] for dimTag in outDimTags_fs]
gmsh.model.addPhysicalGroup(2, tag_fs_surface, tag=MATERIAL[0]['tag'], name=MATERIAL[0]['name'])
print(f"  Tag {MATERIAL[0]['tag']} ('{MATERIAL[0]['name']}'): {tag_fs_surface}")

# Superfície 'pml'
tag_pml_surface = [dimTag[1] for dimTag in outDimTags_pml]
gmsh.model.addPhysicalGroup(2, tag_pml_surface, tag=MATERIAL[1]['tag'], name=MATERIAL[1]['name'])
print(f"  Tag {MATERIAL[1]['tag']} ('{MATERIAL[1]['name']}'): {tag_pml_surface}")

# Definir tamanho e ordem dos elementos
gmsh.option.setNumber("Mesh.MeshSizeMax", h)
gmsh.option.setNumber("Mesh.MeshSizeMin", h)
gmsh.model.mesh.generate(2)
gmsh.model.mesh.setOrder(1)

# Visualizar a malha no ambiente Gmsh (opcional)
gmsh.fltk.run()

# Finalizar o Gmsh
gmsh.finalize()