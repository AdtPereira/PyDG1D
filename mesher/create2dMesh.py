import gmsh


def createSquareK14(filename: str):
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
    gmsh.finalize()