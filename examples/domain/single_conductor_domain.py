import gmsh
import os

# ----------------------------
# Configurações iniciais
# ----------------------------

os.system('cls' if os.name == 'nt' else 'clear')
gmsh.initialize()
gmsh.model.add("final_geometry")
factory = gmsh.model.occ

# ----------------------------
# Parâmetros da geometria
# ----------------------------
x0, y0 = 1.0, 1.0     # semi-lados do retângulo interno
L = 1.0               # camada externa
ra = 0.5              # raio do disco central
h = 0.2               # tamanho de discretização

# ----------------------------
# Geometria base
# ----------------------------
outer = factory.addRectangle(-(x0 + L), -(y0 + L), 0, 2*(x0 + L), 2*(y0 + L))
inner = factory.addRectangle(-x0, -y0, 0, 2*x0, 2*y0)
disk  = factory.addDisk(0, 0, 0, ra, ra)

# ----------------------------
# Operações booleanas (sem sobreposição)
# ----------------------------
region_inner, _ = factory.cut([(2, inner)], [(2, disk)], removeObject=True, removeTool=True)
region_outer, _ = factory.cut([(2, outer)], [(2, inner)], removeObject=True, removeTool=False)
region_core  = [(2, disk)]

# ----------------------------
# Sincronização
# ----------------------------
factory.synchronize()

# ----------------------------
# Malha
# ----------------------------
gmsh.option.setNumber("Mesh.MeshSizeMax", h)
gmsh.option.setNumber("Mesh.MeshSizeMin", h)
gmsh.model.mesh.generate(2)
gmsh.model.mesh.setOrder(1)
gmsh.fltk.run()
gmsh.finalize()
