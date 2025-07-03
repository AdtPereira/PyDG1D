import os
import gmsh
import numpy as np
from pathlib import Path
from scipy.constants import mu_0, epsilon_0

# Limpar o terminal
os.system('cls' if os.name == 'nt' else 'clear')

# --- Parâmetros Físicos ---
FREQ = 3e8
OMEGA = 2 * np.pi * FREQ
K0 = OMEGA * np.sqrt(mu_0 * epsilon_0)
WAVELENGTH = 2 * np.pi / K0

# --- Parâmetros Geométricos e de Malha ---
h = WAVELENGTH / 0.10
L = WAVELENGTH / 2
x0 = WAVELENGTH

print(f"Wavelength: {WAVELENGTH:.4f} m")
print(f"Mesh element size (h): {h:.4f} m")

# Inicializar o Gmsh
gmsh.initialize()
gmsh.model.add("cubeK168_upml")
factory = gmsh.model.occ

# --- Criação da Geometria ---
total_field_domain = factory.addBox(-x0, -x0, -x0, 2*x0, 2*x0, 2*x0)
scattered_field_boundary = factory.addBox(-x0 - L, -x0 - L, -x0 - L, 2*(x0 + L), 2*(x0 + L), 2*(x0 + L))
pml_boundary = factory.addBox(-x0 - 2*L, -x0 - 2*L, -x0 - 2*L, 2*(x0 + 2*L), 2*(x0 + 2*L), 2*(x0 + 2*L))

pml_map, _ = factory.cut([(3, pml_boundary)], [(3, scattered_field_boundary)], removeTool=False)
scattered_field_map, _ = factory.cut([(3, scattered_field_boundary)], [(3, total_field_domain)], removeTool=False)

factory.synchronize()

# --- Definição dos Grupos Físicos ---
# Extrai a tag numérica de cada volume criado.
# O resultado das operações de corte pode ter mais de um volume, mas aqui sabemos que é apenas um.
pml_tag = pml_map[0][1]
scattered_field_tag = scattered_field_map[0][1]
# A tag do domínio total-field é a original, que geralmente é 1
total_field_tag = total_field_domain

# Atribui cada volume a um grupo físico com um nome e uma tag numérica única.
# A dimensão é 3 para volumes.
gmsh.model.addPhysicalGroup(3, [total_field_tag], tag=301, name="Total Field")
gmsh.model.addPhysicalGroup(3, [scattered_field_tag], tag=302, name="Scattered Field")
gmsh.model.addPhysicalGroup(3, [pml_tag], tag=303, name="PML")

# --- Malha ---
gmsh.option.setNumber("Mesh.MeshSizeMax", h)
gmsh.option.setNumber("Mesh.MeshSizeMin", h)
gmsh.model.mesh.generate(3)

# --- Saída de Arquivos ---
output_dir = Path("examples/domain/upml/")
output_dir.mkdir(parents=True, exist_ok=True)
msh_file = output_dir / "cubeK168.msh"
vtk_file = output_dir / "cubeK168.vtk"

gmsh.write(str(msh_file))
gmsh.write(str(vtk_file))
print(f"\nMalha salva em {msh_file}")

# --- Visualização ---
print("Iniciando a GUI do Gmsh...")
gmsh.fltk.run()

# Finalizar o Gmsh
gmsh.finalize()