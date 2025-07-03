import sys
import os

import gmsh
import numpy as np
from pathlib import Path
from scipy.constants import mu_0, epsilon_0

# Limpar o terminal
import os   
os.system('cls' if os.name == 'nt' else 'clear')

# Parâmetros físicos
OMEGA = 2 * np.pi * 3E8
K0 = OMEGA * np.sqrt(mu_0 * epsilon_0)
WAVELENGTH = 2 * np.pi / K0

h = 10*WAVELENGTH      # Tamanho do elemento finito
L = WAVELENGTH/2    # Lado do cubo externo
ra = WAVELENGTH/2   # Raio da esfera
x0 = WAVELENGTH     # Lado do cubo interno

print(f"\nTamanho do elemento: {h:.2e}")
print(f"Lado do cubo externo: {L:.2e} m")
print(f"Raio da esfera: {ra:.2e} m")
print(f"Lado do cubo interno: {x0:.2e} m")

# Inicializar o Gmsh
gmsh.initialize()
gmsh.model.add("cavity_pml")
factory = gmsh.model.occ

# PML cúbicas Inferior Esquerda
region_a = factory.addBox(-x0-L, -x0-L, +x0, L, L, L)
region_b = factory.addBox(-x0-L, -x0-L, -x0-L, L, L, L)

# PML cúbicas Inferior Direita
region_c = factory.addBox(x0, -x0-L, +x0, L, L, L)
region_d = factory.addBox(x0, -x0-L, -x0-L, L, L, L)

# PML cúbicas Superior Direita
region_e = factory.addBox(x0, x0, +x0, L, L, L)
region_f = factory.addBox(x0, x0, -x0-L, L, L, L)

# PML cúbicas Inferior Esquerda
region_g = factory.addBox(-x0-L, x0, +x0, L, L, L)
region_h = factory.addBox(-x0-L, x0, -x0-L, L, L, L)

# PML longitudinais Inferiores
region_1 = factory.addBox(x0, -x0-L, -x0, L, L, 2*x0)
region_2 = factory.addBox(-x0-L, -x0-L, -x0, L, L, 2*x0)

# PML longitudinais Superiores
region_3 = factory.addBox(x0, x0, -x0, L, L, 2*x0)
region_4 = factory.addBox(-x0-L, x0, -x0, L, L, 2*x0)

# PML horizontais Inferiores
region_5 = factory.addBox(-x0, -x0-L, x0, 2*x0, L, L)
region_6 = factory.addBox(-x0, -x0-L, -x0-L, 2*x0, L, L)

# PML horizontais superiores
region_7 = factory.addBox(-x0, x0, x0, 2*x0, L, L)
region_8 = factory.addBox(-x0, x0, -x0-L, 2*x0, L, L)

# PML verticais direita
region_9 = factory.addBox(x0, -x0, x0, L, 2*x0, L)
region_10 = factory.addBox(x0, -x0, -x0-L, L, 2*x0, L)

# PML verticais esquerda
region_11 = factory.addBox(-x0-L, -x0, x0, L, 2*x0, L)
region_12 = factory.addBox(-x0-L, -x0, -x0-L, L, 2*x0, L)

# PML Piso e Teto
region_13 = factory.addBox(-x0, -x0-L, -x0, 2*x0, L, 2*x0)
region_14 = factory.addBox(-x0, x0, -x0, 2*x0, L, 2*x0)

# PML Parede Direita e Esquerda
region_15 = factory.addBox(x0, -x0, -x0, L, 2*x0, 2*x0)
region_16 = factory.addBox(-x0-L, -x0, -x0, L, 2*x0, 2*x0)

# PML Parede Frontal e Traseira
region_17 = factory.addBox(-x0, -x0, x0, 2*x0, 2*x0, L)
region_18 = factory.addBox(-x0, -x0, -x0-L, 2*x0, 2*x0, L)

# Criar região do espaço livre (cubo interno)
region_fs = factory.addBox(-x0, -x0, -x0, 2*x0, 2*x0, 2*x0)

# Criar região do espalhador (esfera)
sphere = factory.addSphere(0, 0, 0, ra)

# Regiões PML em termos de coordenadas
PML_xyz = [region_a, region_b, region_c, region_d, region_e, region_f, region_g, region_h]
PML_xy = [region_1, region_2, region_3, region_4]
PML_yz = [region_5, region_6, region_7, region_8]
PML_xz = [region_9, region_10, region_11, region_12]
PML_y = [region_13, region_14]
PML_x = [region_15, region_16]
PML_z = [region_17, region_18]

# Fragmentar todas as regiões para garantir interfaces conformais
PML_list = PML_xy + PML_yz + PML_xz + PML_y + PML_x + PML_z + PML_xyz
objectDimTags = [(3, region_fs)] + [(3, item) for item in PML_list]
outDimTags, _ = factory.fragment(objectDimTags, objectDimTags)

# Subtrair a esfera do espaço livre
outDimTags_omega_s, _ = factory.cut([(3, region_fs)], [(3, sphere)], removeTool=True)

# Sincronizar
factory.synchronize()

# Definir ordem dos elementos
gmsh.option.setNumber("Mesh.MeshSizeMax", h)
gmsh.option.setNumber("Mesh.MeshSizeMin", h)
gmsh.model.mesh.generate(3)

# Adicionar o diretório raiz do projeto ao sys.path
project_root = Path().resolve().parent  
sys.path.append(str(project_root))
os.makedirs("examples/domain/upml/", exist_ok=True)
gmsh.write("examples/domain/upml/cubeK494.msh")
gmsh.write("examples/domain/upml/cubeK494.vtk")

# Visualizar a malha
gmsh.fltk.run()

# Finalizar o Gmsh
gmsh.finalize()
