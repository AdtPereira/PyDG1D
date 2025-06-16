import os
import gmsh

# Limpar o terminal
os.system('cls' if os.name == 'nt' else 'clear')

WAVELENGTH = 1.0  # Comprimento de onda em metros

BOUNDARY = [{'tag': 101, 'type': 'Dirichlet', 'value': None, 'name': 'circular_scatterer'}]

MATERIAL = [{'tag': 201, 'name': 'free_space', 'relative_magnetic_permeability': 1, 'relative_electric_permittivity': 1},
            {'tag': 301, 'name': 'PML_a','relative_magnetic_permeability': 1, 'relative_electric_permittivity': 1},
            {'tag': 302, 'name': 'PML_b','relative_magnetic_permeability': 1, 'relative_electric_permittivity': 1},
            {'tag': 303, 'name': 'PML_c','relative_magnetic_permeability': 1, 'relative_electric_permittivity': 1},
            {'tag': 304, 'name': 'PML_d','relative_magnetic_permeability': 1, 'relative_electric_permittivity': 1},
            {'tag': 401, 'name': 'PML_I','relative_magnetic_permeability': 1, 'relative_electric_permittivity': 1},
            {'tag': 402, 'name': 'PML_II','relative_magnetic_permeability': 1, 'relative_electric_permittivity': 1},
            {'tag': 403, 'name': 'PML_III','relative_magnetic_permeability': 1, 'relative_electric_permittivity': 1},
            {'tag': 404, 'name': 'PML_IV','relative_magnetic_permeability': 1, 'relative_electric_permittivity': 1}]

GEOMETRY = {'h': WAVELENGTH*0.5, 'L': WAVELENGTH, 'ra': WAVELENGTH*0.0, 'x0': 1.0*WAVELENGTH}

# Inicializar o Gmsh
gmsh.initialize()
gmsh.model.add("rectangular_pml")
factory = gmsh.model.occ

# Dimensões do domínio 
h = GEOMETRY['h']       # Tamanho do elemento finito
L = GEOMETRY['L']       # Espessura da PML
ra = GEOMETRY['ra']     # Raio do espalhador
x0 = GEOMETRY['x0']     # Fronteira interna da PML
xa, xb = x0, x0 + L     # Fronteiras externas da PML

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

# Criar região do espalhador, omega_s
# disk = factory.addDisk(0, 0, 0, ra, ra)

# Subtrair omega_s de omega_fs
# outDimTags_omega_s, _ = factory.cut([(2, region_fs)], [(2, disk)], removeTool=True)
    
# Sincronizar após o corte do retângulo interno
factory.synchronize()

# Obter o contorno (curva, dim=1) de gamma_s
# DimTags_free_space = gmsh.model.getBoundary(outDimTags_omega_s, oriented=True, recursive=False)

# TAGs de gamma_s
# TagList_scatterer = [-tag[1] for tag in DimTags_free_space if tag[1] < 0]

# TAGs de gamma_fs
# tagList_fs = [tag[1] for tag in DimTags_free_space if tag[1] > 0]

# Adicionar grupos físicos para curvas (Dim=1)
# for i, CurveTagList in enumerate([TagList_scatterer]):
#     gmsh.model.addPhysicalGroup(1, CurveTagList, tag=BOUNDARY[i]['tag'], name=BOUNDARY[i]['name'])

# Adicionar grupos físicos para superfícies (Dim=2)	    
TagList_surfaces = [region_fs, region_a, region_b, region_c, region_d, region_i, region_ii, region_iii, region_iv]
for i, SurfaceList in enumerate(TagList_surfaces):
    gmsh.model.addPhysicalGroup(2, [SurfaceList], tag=MATERIAL[i]['tag'], name=MATERIAL[i]['name'])

# Definir ordem dos elementos
gmsh.option.setNumber("Mesh.MeshSizeMax", h)
gmsh.option.setNumber("Mesh.MeshSizeMin", h)
gmsh.model.mesh.generate(2)
gmsh.model.mesh.setOrder(1)

gmsh.fltk.run()
gmsh.finalize()