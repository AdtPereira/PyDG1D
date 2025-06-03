import gmsh
import os

# Limpar o terminal
os.system('cls' if os.name == 'nt' else 'clear')

# Define os parâmetros de entrada
FINITE_ELEMENT = ("Triangle", 1)
BOUNDARY = [{'tag': 101,
             'type': 'Neumann',
             'value': 0.0,
             'name': 'Ez_0'}]
MATERIAL = [{'tag': 201,
             'name': 'free_space',
             'relative_magnetic_permeability': 1,
             'relative_electric_permittivity': 1}]

element, order = FINITE_ELEMENT
h = 1

# Dimensões do guia de onda retangular
a, b = 2, 2

# Inicializar o Gmsh
gmsh.initialize()
gmsh.model.add("rectangular_guide")

# Criar superfície retangular
TagSurface = gmsh.model.occ.addRectangle(-a/2, -b/2, 0, a, b)
gmsh.model.occ.synchronize()
gmsh.option.setNumber("Mesh.MeshSizeMin", h)
gmsh.model.mesh.generate(dim=2)
gmsh.model.mesh.setOrder(order)

# Obter os contornos (curvas, dim=1) de cada superfície
outDimTags = gmsh.model.getBoundary([(2, TagSurface)], oriented=True, recursive=False)

# Exibir os TAGs das curvas associadas a cada contorno
tagList_boundary = [Dimtags[1] for Dimtags in outDimTags]

# Definindo as curvas de contorno de Dirichlet (dim=1)
gmsh.model.addPhysicalGroup(dim=1, tags=tagList_boundary, tag=BOUNDARY[0]['tag'], name=BOUNDARY[0]['name'])

# Adicionar grupos físicos para Dim=2 (superfícies)
gmsh.model.addPhysicalGroup(dim=2, tags=[TagSurface], tag=MATERIAL[0]['tag'], name=MATERIAL[0]['name'])

# Visualizar a malha no ambiente Gmsh (opcional)
gmsh.fltk.run()

# Finalizar o Gmsh
gmsh.finalize()