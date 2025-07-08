import os
import gmsh

if __name__ == "__main__":
    # Clear terminal
    os.system('cls' if os.name == 'nt' else 'clear')

    # --- Parâmetros Geométricos e de Malha ---
    DIM = 3                         # Dimensão do problema
    L = 1.0                         # Largura da camada da PML
    h = 4.0                         # Tamanho máximo do elemento da malha
    Lx = 2.0                        # Dimensão total do domínio na direção x
    x0 = 1.0                        # semi-lados do retângulo interno - Domínio Físico

    # 1. Inicializa o GMSH
    gmsh.initialize()

    gmsh.model.add("cubeK96_upml")
    factory = gmsh.model.occ

    # --- Criação da Geometria ---
    tfz = factory.addBox(-x0, -x0, -x0, 2*x0, 2*x0, 2*x0)
    pml = factory.addBox(-Lx, -Lx, -Lx, 2*Lx, 2*Lx, 2*Lx)
    pml_domain, _ = factory.cut([(3, pml)], [(3, tfz)], removeTool=False)
    factory.synchronize()

    # --- Definição dos Grupos Físicos ---
    pmlTag = pml_domain[0][1]
    tfzTag = tfz

    # Atribui cada volume a um grupo físico com um nome e uma tag numérica única.
    gmsh.model.addPhysicalGroup(3, [tfzTag], tag=301, name="TFZ")
    gmsh.model.addPhysicalGroup(3, [pmlTag], tag=303, name="PML")

    # --- Malha ---
    gmsh.option.setNumber("Mesh.MeshSizeMax", h)
    gmsh.option.setNumber("Mesh.MeshSizeMin", h)
    gmsh.model.mesh.generate(3)

    # 6. Exibe a malha e Finaliza o GMSH
    gmsh.fltk.run()
    gmsh.finalize()
