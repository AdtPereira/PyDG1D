import gmsh


def SquareK14Domain(filename: str):
    gmsh.initialize()
    gmsh.model.add("SquareK14")

    # 1. Criar superfície retangular
    TagSurface = gmsh.model.occ.addRectangle(-1, -1, 0, 2, 2)
    gmsh.model.occ.synchronize()
    gmsh.option.setNumber("Mesh.MeshSizeMin", 1.0)    
    gmsh.model.mesh.generate(dim=2)
    gmsh.model.mesh.setOrder(1)

    # Obter os contornos (curvas, dim=1) de cada superfície
    outDimTags = gmsh.model.getBoundary([(2, TagSurface)], oriented=True, recursive=False)
    tagList_boundary = [Dimtags[1] for Dimtags in outDimTags]

    # Adicionar grupos físicos para Dim=1 (curvas)
    gmsh.model.addPhysicalGroup(dim=1, tags=tagList_boundary, tag=101, name='outerBoundary')

    # Adicionar grupos físicos para Dim=2 (superfícies)
    gmsh.model.addPhysicalGroup(dim=2, tags=[TagSurface], tag=201, name='freeSpace')

    gmsh.write(filename)
    # gmsh.fltk.run()
    # gmsh.finalize()


def SquareK16Domain(problem, filename: str):
    # Problem dimensions
    dim = problem['DIM']            # Dimensionality of the problem
    h = problem['domain']['h']      # Mesh size
    Lx = problem['domain']['Lx']    # Half-length in x-direction
    Ly = problem['domain']['Ly']    # Half-length in y-direction
    L = problem['pml']['L']         # PML layer width
    x0 = Lx - L                     # semi-lados do retângulo interno - Domínio Físico

    gmsh.initialize()
    gmsh.model.add("SquareK16")
    factory = gmsh.model.occ

    # 1. Criar superfície retangular
    # --- Criação da Geometria ---
    tfz = factory.addRectangle(-x0, -x0, 0, 2*x0, 2*x0)
    pml = factory.addRectangle(-Lx, -Ly, 0, 2*Lx, 2*Ly)
    pml_domain, _ = factory.cut([(2, pml)], [(2, tfz)], removeTool=False)
    gmsh.model.occ.synchronize()
    
    # --- Definição dos Grupos Físicos ---
    # Extrai a tag numérica de cada volume criado.
    # O resultado das operações de corte pode ter mais de um volume, mas aqui sabemos que é apenas um.
    pmlTag = pml_domain[0][1]
    # A tag do domínio total-field é a original, que geralmente é 1
    tfzTag = tfz

    # Obter os contornos (curvas, dim=1) de cada superfície
    outDimTags = gmsh.model.getBoundary([(2, pmlTag)], oriented=True, recursive=False)
    tfzInterfaceTags = [abs(Dimtags[1]) for Dimtags in outDimTags if Dimtags[1] < 0]
    pmlInterfaceTags = [Dimtags[1] for Dimtags in outDimTags if Dimtags[1] > 0]

    # Adicionar grupos físicos para Dim=1 (curvas)
    gmsh.model.addPhysicalGroup(dim=1, tags=pmlInterfaceTags, tag=101, name='outerBoundary')
    gmsh.model.addPhysicalGroup(dim=1, tags=tfzInterfaceTags, tag=102, name='tfzInterface')

    # Adicionar grupos físicos para Dim=2 (superfícies)
    gmsh.model.addPhysicalGroup(dim=2, tags=[pmlTag], tag=201, name='pml')
    gmsh.model.addPhysicalGroup(dim=2, tags=[tfzTag], tag=202, name='tfz')
    
    # --- Malha ---
    gmsh.option.setNumber("Mesh.MeshSizeMax", h)
    gmsh.option.setNumber("Mesh.MeshSizeMin", h)
    gmsh.model.mesh.generate(dim)
    gmsh.model.mesh.setOrder(1)

    # gmsh.fltk.run()
    gmsh.write(filename)