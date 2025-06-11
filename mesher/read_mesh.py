# %%
import gmsh


# %% [markdown]
# # Gmsh application programming interface (API)
# 
# Gmsh 4.13.1: A three-dimensional finite element mesh generator with built-in pre- and post-processing facilities. https://gmsh.info/
# 
# The Gmsh application programming interface (API) allows to integrate the Gmsh library in external applications written in C++, C, Python, Julia or Fortran.
# 
# References:  
# 1. Geometry basics, elementary entities, physical groups. Tutorial t1. https://gmsh.info/doc/texinfo/gmsh.html#t1
# 2. https://gitlab.onelab.info/gmsh/gmsh/blob/gmsh_4_13_1/tutorials/python/x1.py#L33
# 3. https://gitlab.onelab.info/gmsh/gmsh/-/tree/gmsh_4_13_1/examples/api
# 4. https://gitlab.onelab.info/gmsh/gmsh/blob/gmsh_4_13_1/api/gmsh.py
# 5. https://gmsh.info/doc/texinfo/gmsh.html#Namespace-gmsh_002fmodel_002fgeo
# 
# ## Namespace ``gmsh/model``: model functions
# 
# Geometrical data is made of elementary model 'entities', called 'points' (entities of dimension 0), `curves' (entities of dimension 1), 'surfaces' (entities of dimension 2) and 'volumes' (entities of dimension 3). Elementary model entities are identified by their dimension and by a 'tag': a strictly positive identification number. 'Physical groups' are collections of model entities and are also identified by their dimension and by a tag.
# 
# ### ``gmsh/model/getPhysicalGroups``
# 
# > Get all the physical groups in the current model. If ``dim`` is >= 0, return only the entities of the specified dimension (e.g. physical points if ``dim`` == 0). **The entities are returned as a vector of (``dim``, ``tag``) pairs.**
# 
# > Input: dim = -1 (integer)  
# > Output: dimTags (vector of pairs of integers)  
# > Return: -
# 
# ### ``gmsh/model/getEntitiesForPhysicalGroup``
# 
# > Get the tags of the model entities making up the physical group of dimension ``dim`` and tag ``tag``.
# 
# > Input: dim (integer), tag (integer)  
# > Output: tags (vector of integers)  
# > Return: -
# 
# ### ``gmsh/model/getEntities``
# 
# > Get all the entities in the current model. A model entity is represented by two integers: its dimension (dim == 0, 1, 2 or 3) and its tag (its unique, strictly positive identifier). If ``dim`` is >= 0, return only the entities of the specified dimension (e.g. points if ``dim`` == 0). The entities are returned as a vector of (dim, tag) pairs.
# 
# > Input: ``dim`` = -1 (integer)  
# > Output: ``dimTags`` (vector of pairs of integers)  
# > Return: -
# 
# ## Namespace gmsh/model/mesh: mesh functions
# 
# ### ``gmsh/model/mesh/getNodes``
# 
# > Get the nodes classified on the entity of dimension ``dim`` and tag ``tag``. If ``tag`` < 0, get the nodes for all entities of dimension ``dim``. If ``dim`` and ``tag`` are _negative_, get all the nodes in the mesh. ``nodeTags`` contains the node tags (their unique, strictly positive identification numbers). coord is a vector of length 3 times the length of nodeTags that contains the x, y, z coordinates of the nodes, concatenated: [n1x, n1y, n1z, n2x, ...]. If dim >= 0 and return ParamtricCoord is set, parametricCoord contains the parametric coordinates ([u1, u2, ...] or [u1, v1, u2, ...]) of the nodes, if available. The length of parametricCoord can be 0 or dim times the length of nodeTags. If includeBoundary is set, also return the nodes classified on the boundary of the entity (which will be reparametrized on the entity if dim >= 0 in order to compute their parametric coordinates).
# 
# > Input: dim = -1 (integer), tag = -1 (integer), includeBoundary = False (boolean), returnParametricCoord = True (boolean)  
# > Output: nodeTags (vector of sizes), coord (vector of doubles), parametricCoord (vector of doubles)  
# > Return: -

# %%
def get_node_coordinates():
    # Obter as coordenadas dos nós
    nodeTags, nodeCoords, nodeParams = gmsh.model.mesh.getNodes()
    nodes = [(nodeCoords[i], nodeCoords[i + 1]) for i in range(0, len(nodeCoords), 3)]
    return nodes

# %% [markdown]
# ### ``gmsh/model/mesh/getElements``
# 
# > Get the elements classified on the entity of dimension `dim` and tag ``tag``. If ``tag`` < 0, get the elements for all entities of dimension ``dim``. If ``dim`` and ``tag`` are negative, get all the elements in the mesh. ``elementTypes`` contains the MSH types of the elements (e.g. ``2`` for 3-node triangles: see ``getElementProperties`` to obtain the properties for a given element type). ``elementTags`` is a vector of the same length as ``elementTypes``; each entry is a vector containing the tags (unique, strictly positive identifiers) of the elements of the corresponding type. ``nodeTags`` is also a vector of the same length as ``elementTypes``; each entry is a vector of length equal to the number of elements of the given type times the number N of nodes for this type of element, that contains the node tags of all the elements of the given type, concatenated: [e1n1, e1n2, ..., e1nN, e2n1, ...].
# 
# > Input: ``dim`` = -1 (integer), ``tag`` = -1 (integer)
# > Output: ``elementTypes`` (vector of integers), ``elementTags`` (vector of vectors of sizes), ``nodeTags`` (vector of vectors of sizes)
# > Return: -
# 
# ### `gmsh/model/mesh/getElementProperties`
# 
# > Get the properties of an element of type `elementType`: its name (``elementName``), dimension (``dim``), order (``order``), number of nodes (``numNodes``), local coordinates of the nodes in the reference element (``localNodeCoord`` vector, of length ``dim`` times ``numNodes``) and number of primary (first order) nodes (``numPrimaryNodes``).
# 
# > Input: ``elementType`` (integer)  
# > Output: ``elementName`` (string), ``dim`` (integer), ``order`` (integer), ``numNodes`` (integer), ``localNodeCoord`` (vector of doubles), ``numPrimaryNodes`` (integer)  
# > Return: -

# %% [markdown]
# # `get_conn()`

# %%
def get_conn():
    # Obter os elementos da malha
    elemTypes, elemTags, elemNodeTags = gmsh.model.mesh.getElements(dim=2)
    conn = []

    for elemType, elemNode in zip(elemTypes, elemNodeTags):
        # Obter as propriedades do elemento
        name, _, _, nodes_per_element, _, _ = gmsh.model.mesh.getElementProperties(elemType)
        numElements = len(elemNode) // nodes_per_element
        for i in range(numElements):
            conn_i = elemNode[nodes_per_element * i: nodes_per_element * (i + 1)]
            conn.append(conn_i)

    # Retornar a matriz de conectividade com listas para cada elemento
    return [element.tolist() for element in conn]

# %% [markdown]
# # `get_cell_data()`

# %%
def get_cell_data(MATERIAL, dim=2):

    # 1. Criar um mapa entre cada entidade física (grupo físico) e os elementos correspondentes
    physical_groups_map = {}
    for material in MATERIAL:
        geo_entities = gmsh.model.getEntitiesForPhysicalGroup(dim, material['tag'])
        for EntityTag in geo_entities:
            mesh_elements = gmsh.model.mesh.getElements(dim, EntityTag)
            physical_groups_map.update({element: material for element in mesh_elements[1][0]})

    # 2. Obter os elementos da malha
    cell_data = {}
    elemTypes, elemTags, elemNodeTags = gmsh.model.mesh.getElements(dim)
    
    for elemType, elemTag, elemNode in zip(elemTypes, elemTags, elemNodeTags):        
        # Obter as propriedades do elemento
        _, _, _, nodes_per_element, _, _ = gmsh.model.mesh.getElementProperties(elemType)
        
        # Número de elementos
        Ne = len(elemNode) // nodes_per_element

        for i in range(Ne):
            # Obter a conectividade do elemento
            conn_node = elemNode[nodes_per_element * i: nodes_per_element * (i + 1)].tolist()

            # Adicionar informações do elemento ao dicionário
            cell_data[i + 1] = {
                'conn': conn_node,
                'conn_sorted': sorted(conn_node),
                'conn_edge': None,
                'geo': {'centroid': None, 'dim': None},
                'contour': {'type': None, 'conn_contour': None},
                'material': physical_groups_map.get(elemTag[i], None)
            }

    return cell_data

# %% [markdown]
# # `get_new_cell_data()`

# %%
def get_new_cell_data(MATERIAL, problem_dim):

    # 1. Criar um mapa de arestas
    gmsh.model.mesh.createEdges()
    
    # Obter as arestas da malha
    edgeTags, edgeCoords = gmsh.model.mesh.getAllEdges()

    # Criar edge_mapping diretamente sem edge_mapping
    edge_mapping = {tuple(sorted([edgeCoords[2*i], edgeCoords[2*i+1]])): tag
                        for i, tag in enumerate(sorted(edgeTags))}

    # 2. Criar o dicionário mesh_data['cell]
    cell_data = {}
    for material in MATERIAL:
        # Obter as entidades físicas (grupo físico) associadas ao material
        MaterialEntitiesTags = gmsh.model.getEntitiesForPhysicalGroup(problem_dim, tag=material['tag'])
        
        for EntityTag in MaterialEntitiesTags:
            # Obter os elementos da malha
            elemTypes, elemTags, elemNodes = gmsh.model.mesh.getElements(problem_dim, tag=EntityTag)
            
            for elemType, elemTag, elemNode in zip(elemTypes, elemTags, elemNodes):            
                # Obter as propriedades do elemento
                _, _, _, nodes_per_element, _, _ = gmsh.model.mesh.getElementProperties(elemType)

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
                    
                    # Adicionar ao dicionário
                    cell_data[Tag] = {
                        'tag': Tag,
                        'conn': conn_node,
                        'conn_sorted': conn_std,
                        'conn_edge': conn_edge,
                        'conn_face': conn_face,
                        'geo': {'centroid': None, 'dim': None},
                        'contour': {'type': None, 'conn_contour': None},
                        'material': material}
                    
    # 3. Reordenar chave de cell_data
    cell_data = {i + 1: cell_data[Tag] for i, Tag in enumerate(cell_data)}
                    
    return cell_data

# %% [markdown]
# # `get_boundary_data()`

# %%
def get_boundary_data(BOUNDARY, problem_dim):

    # 1. Criar um mapa de arestas
    gmsh.model.mesh.createEdges()
    
    # Obter as arestas da malha
    edgeTags, edgeNodes = gmsh.model.mesh.getAllEdges()

    # Criar edge_mapping diretamente sem edge_mapping
    edge_mapping = {tuple(sorted([edgeNodes[2*i], edgeNodes[2*i+1]])): tag
                        for i, tag in enumerate(sorted(edgeTags))}

    # 3. Criar o dicionário mesh_data['boundary']
    boundary_data = {}
    for bc in BOUNDARY:
        # Obter as entidades físicas (grupo físico) associadas ao contorno
        BoundaryEntitiesTags = gmsh.model.getEntitiesForPhysicalGroup(problem_dim-1, tag=bc['tag'])

        for EntityTag in BoundaryEntitiesTags:
            # Obter os elementos da malha
            elemTypes, elemTags, elemNodes = gmsh.model.mesh.getElements(problem_dim-1, tag=EntityTag)
            
            for elemType, elemTag, elemNode in zip(elemTypes, elemTags, elemNodes):            
                # Obter as propriedades do elemento
                _, _, _, nodes_per_element, _, _ = gmsh.model.mesh.getElementProperties(elemType)

                # Criar dicionário associando cada elemento à sua conectividade de nós
                for i, Tag in enumerate(elemTag):                    
                    # Obter a conectividade do elemento
                    conn_node = elemNode[nodes_per_element * i: nodes_per_element * (i + 1)].tolist()
                    conn_std = sorted(conn_node)

                    # Triangular element
                    if nodes_per_element == 3:
                        conn_edge = [
                            edge_mapping[(conn_std[0], conn_std[1])],  # e1: 1 -> 2
                            edge_mapping[(conn_std[0], conn_std[2])],  # e2: 1 -> 3
                            edge_mapping[(conn_std[1], conn_std[2])]   # e3: 2 -> 3
                        ]

                    # Adicionar ao dicionário
                    boundary_data[Tag] = {
                        'tag': Tag,
                        'conn': conn_node,
                        'conn_sorted': conn_std,
                        'conn_edge': conn_edge,
                        'geo': {'centroid': None, 'dim': None},
                        'contour': {'type': None, 'conn_contour': None},
                        'boundary': bc}
    
    # 3. Reordenar chave de boundary_data
    boundary_data = {i + 1: boundary_data[Tag] for i, Tag in enumerate(boundary_data)}

    return boundary_data

# %% [markdown]
# # `get_nodes_data()`

# %%
def get_nodes_data(BOUNDARY, problem_dim):
    # 1. Dicionário para Mapeamento inicial de nós:
    # Todos os nós começam com a condição de contorno "Free" e valor None.
    NodeTags, NodeCoords, _ = gmsh.model.mesh.getNodes()
    node_mapping = {tag: {'tag': None, 'type': 'Free', 'value': None, 'name': 'free_node'}
                    for tag in NodeTags}

    # 2. Atualização do mapeamento de nós:
    # Para cada condição de contorno em BOUNDARY, os nós associados ao grupo físico correspondente 
    # são atualizados com o tipo (type) e valor (value) dessa condição.
    for bc in BOUNDARY:
        # Obtenha os nós associados ao grupo físico especificado no bc
        bc_NodeTags, _ = gmsh.model.mesh.getNodesForPhysicalGroup(problem_dim-1, tag=bc['tag'])
        
        # Atualiza o mapeamento de nós com a condição de contorno correspondente
        for node in bc_NodeTags:
            node_mapping[node] = {
                'tag': bc['tag'], 'type': bc['type'], 'value': bc['value'], 'name': bc['name']}

    # 3. Estrutura final:
    # A lista dict_nodes contém informações completas sobre cada nó, incluindo suas coordenadas 
    # globais e as condições de contorno associadas.
    nodes_data = {
        node: {"xg": (NodeCoords[3*i], NodeCoords[3*i + 1], NodeCoords[3*i + 2]),
                "bc": node_mapping[node]} 
        for i, node in enumerate(NodeTags)}

    return nodes_data

# %% [markdown]
# # `get_new_nodes_data()`

# %%
def get_new_nodes_data(BOUNDARY, INTERFACES, dim=1):
    # 1. Dicionário para Mapeamento inicial de nós:
    # Todos os nós começam com a condição de contorno "Free" e valor None.
    NodeTags, NodeCoords, _ = gmsh.model.mesh.getNodes()
    node_bc_map = {tag: {'tag': None, 'type': 'Free', 'value': None, 'name': 'free_node'}
                    for tag in NodeTags}

    # 2. Atualização do mapeamento de nós:
    # Para cada condição de contorno em BOUNDARY, os nós associados ao grupo físico correspondente 
    # são atualizados com o tipo (type) e valor (value) dessa condição.
    for bc in BOUNDARY:
        # Obtenha os nós associados ao grupo físico especificado no bc
        condition_NodeTags, _ = gmsh.model.mesh.getNodesForPhysicalGroup(dim, tag=bc['tag'])
        
        # Atualiza o mapeamento de nós com a condição de contorno correspondente
        for node in condition_NodeTags:
            node_bc_map[node] = {'tag': bc['tag'], 'type': bc['type'], 'value': bc['value'], 'name': bc['name']}

    # Para cada condição de contorno em INTERFACES, os nós associados ao grupo físico correspondente 
    # são atualizados com o tipo (type) e valor (value) dessa condição.
    for bc in INTERFACES:
        # Obtenha os nós associados ao grupo físico especificado no bc
        condition_NodeTags, _ = gmsh.model.mesh.getNodesForPhysicalGroup(dim, tag=bc['tag'])
        
        # Atualiza o mapeamento de nós com a condição de contorno correspondente
        for node in condition_NodeTags:
            node_bc_map[node] = {'tag': bc['tag'], 'type': bc['type'], 'value': bc['value'], 'name': bc['name']}

    # 3. Estrutura final:
    # A lista dict_nodes contém informações completas sobre cada nó, incluindo suas coordenadas 
    # globais e as condições de contorno associadas.
    dict_nodes = {
        node: {
            "xg": (NodeCoords[3*i], NodeCoords[3*i + 1], NodeCoords[3*i + 2]),
            "bc": node_bc_map[node]} 
        for i, node in enumerate(NodeTags)}

    return dict_nodes

# %% [markdown]
# # `get_edge_data()`

# %%
def get_edge_data():
    # Edges and faces are returned for each element as a list of nodes corresponding
    # to the canonical orientation of the edges and faces for a given element type.

    # Gmsh can also identify unique edges and faces (a single edge or face whatever
    # the ordering of their nodes) and assign them a unique tag. This identification
    # can be done internally by Gmsh (e.g. when generating keys for basis
    # functions), or requested explicitly as follows:
    gmsh.model.mesh.createEdges()

    # If all you need is the list of all edges or faces in terms of their nodes, you
    # can also directly call:
    edgeTags, edgeNodes = gmsh.model.mesh.getAllEdges()

    # Create connection between edges and elements
    edge_mapping = {}
    for i, tag in enumerate(edgeTags):
        # Adiciona ao mapeamento usando o TAG da aresta como chave
        # Percorre os nós em pares e associa aos AlledgeTags
        edge_mapping[tag] = sorted([edgeNodes[2*i], edgeNodes[2*i + 1]])

    # 1. Dicionário para Mapeamento inicial de arestas:
    # Todos as arestas começam com a condição de contorno "Free" e valor None.
    
    # A lista dict_nodes contém informações completas sobre cada nó, incluindo suas coordenadas 
    # globais e as condições de contorno associadas.
    edges_data = {
        edge: {
            "conn": edge_mapping[edge],
            "bc": {'tag': None, 'type': 'Free', 'value': None, 'name': 'free_edge'}
            } for edge in edgeTags
    }
    
    # Ordena o dicionário edge_mapping com base nas chaves
    sorted_dict_edges = {key: edges_data[key] for key in sorted(edges_data.keys())}
    
    return sorted_dict_edges

# %% [markdown]
# # `get_new_edge_data()`

# %%
def get_new_edge_data(BOUNDARY, problem_dim):
    
    # 1. Obter as arestas da malha
    gmsh.model.mesh.createEdges()
    edgeTags, edgeNodes = gmsh.model.mesh.getAllEdges()

    # Criar edge_mapping diretamente sem edge_mapping
    edge_mapping = {tuple(sorted([edgeNodes[2*i], edgeNodes[2*i+1]])): tag
                        for i, tag in enumerate(sorted(edgeTags))}
    
    # 2. Dicionário para Mapeamento inicial de arestas:
    # A lista edge_data contém informações completas sobre cada aresta, incluindo suas coordenadas 
    # globais e as condições de contorno associadas.
    edge_data = { 
        edge: {'conn': sorted([edgeNodes[2*i], edgeNodes[2*i + 1]]),
                'bc': {'tag': None, 'type': 'Free', 'value': None, 'name': 'free_edge'}}
                for i, edge in enumerate(sorted(edgeTags))}
    
    # 3. Atualização do mapeamento de arestas:
    # Para cada condição de contorno em BOUNDARY, as arestas associadas ao grupo físico correspondente 
    # são atualizados com o tipo (type) e valor (value) dessa condição.
    for bc in BOUNDARY:
        # Obter as entidades físicas (grupo físico) associadas ao contorno
        BoundaryEntitiesTags = gmsh.model.getEntitiesForPhysicalGroup(problem_dim-1, tag=bc['tag'])
        
        for EntityTag in BoundaryEntitiesTags:
            # Obter os elementos da malha
            elemTypes, elemTags, elemNodes = gmsh.model.mesh.getElements(problem_dim-1, tag=EntityTag)
            
            for elemType, elemTag, elemNode in zip(elemTypes, elemTags, elemNodes):            
                # Obter as propriedades do elemento
                _, _, _, nodes_per_element, _, _ = gmsh.model.mesh.getElementProperties(elemType)

                # Criar dicionário associando cada elemento à sua conectividade de nós
                for i, Tag in enumerate(elemTag):                    
                    # Obter a conectividade ordenada do elemento
                    conn_std = sorted(elemNode[nodes_per_element * i: nodes_per_element * (i + 1)].tolist())

                    # Triangular element
                    if nodes_per_element == 3:
                        conn_edge = [
                            edge_mapping[(conn_std[0], conn_std[1])],  # e1: 1 -> 2
                            edge_mapping[(conn_std[0], conn_std[2])],  # e2: 1 -> 3
                            edge_mapping[(conn_std[1], conn_std[2])]   # e3: 2 -> 3
                        ]
                    
                    # Localizar as arestas correspondentes no mapeamento de arestas
                    for edge in conn_edge:
                        edge_data[edge]['bc'] = {'tag': bc['tag'], 'type': bc['type'], 'value': bc['value'], 'name': bc['name']}    

    # 5. Reordenar as chaves de edge_data
    edge_data = {i + 1: edge_data[Tag] for i, Tag in enumerate(edge_data)}

    return edge_data

# %% [markdown]
# ### ``gmsh/model/getBoundary``
# 
# > Get the boundary of the model entities dimTags, given as a vector of (dim, tag) pairs. Return in outDimTags the boundary of the individual entities (if combined is false) or the boundary of the combined geometrical shape formed by all input entities (if combined is true). Return tags multiplied by the sign of the boundary entity if oriented is true. Apply the boundary operator recursively down to dimension 0 (i.e. to points) if recursive is true.
# 
# > Input: dimTags (vector of pairs of integers), combined = True (boolean), oriented = True (boolean), recursive = False (boolean)
# > Output: outDimTags (vector of pairs of integers)
# > Return: -

# %% [markdown]
# ## `get_boundary_nodes()` 

# %%
def get_boundary_nodes():
    # Obter os grupos físicos
    physical_groups = gmsh.model.getPhysicalGroups()
    boundary_nodes = set()

    # Iterar sobre todos os grupos físicos
    for dim, tag in physical_groups:
        # Verifica se o grupo físico é de dimensão 1 (linhas/arestas) com a tag de 'boundary'
        if dim == 1 and tag == 101:
            entities = gmsh.model.getEntitiesForPhysicalGroup(dim, tag)
            for entity in entities:
                # Obter os nós associados à entidade (aresta)
                node_tags, _, _ = gmsh.model.mesh.getNodes(dim=1, tag=entity)
                boundary_nodes.update(node_tags)  # Adiciona os nós das arestas

                # Obter os pontos (vértices) de cada aresta
                start_node, end_node = gmsh.model.getBoundary([(dim, entity)], oriented=False)
                for vertex in [start_node, end_node]:
                    vertex_nodes, _, _ = gmsh.model.mesh.getNodes(dim=0, tag=vertex[1])
                    boundary_nodes.update(vertex_nodes)  # Adiciona os nós dos vértices

    return sorted(boundary_nodes)

# %% [markdown]
# # `basic_info()`

# %%
def basic_info(dim=2):

    # ---------------------------------------------------------------------------------
    #  Reduced version of the Gmsh Python extended tutorial 1
    #  https://gitlab.onelab.info/gmsh/gmsh/blob/gmsh_4_13_1/tutorials/python/x1.py#L33
    # ---------------------------------------------------------------------------------

    # Print the model name and dimension:
    print('Model ' + gmsh.model.getCurrent() + ' (' +
        str(gmsh.model.getDimension()) + 'D)')

    # Get the mesh nodes for the entity (dim, tag):
    nodeTags, nodeCoords, nodeParams = gmsh.model.mesh.getNodes()

    # Get the mesh elements for the entity (dim, tag):
    elemTypes, elemTags, elemNodeTags = gmsh.model.mesh.getElements(dim)
    
    # Calculate total number of elements
    total_elements = sum(len(tags) for tags in elemTags)


    # List of all edges in terms of their nodes
    gmsh.model.mesh.createEdges()
    edgeTags, edgeNodes = gmsh.model.mesh.getAllEdges()

    # Print Basic Info Data
    print("Info     : %d geometric entities" % len(gmsh.model.getEntities()))
    print("Info     : %d physical groups" % len(gmsh.model.getPhysicalGroups()))
    print("Info     : %d nodes in total" % len(nodeTags))
    print("Info     : %d edges in total" % len(edgeTags))
    print(f"Info     : %d {dim}-D elements in total" % total_elements)

# %% [markdown]
# # `complete_info()`

# %%
def complete_info():

    # ---------------------------------------------------------------------------------
    #  Reduced version of the Gmsh Python extended tutorial 1
    #  https://gitlab.onelab.info/gmsh/gmsh/blob/gmsh_4_13_1/tutorials/python/x1.py#L33
    # ---------------------------------------------------------------------------------

    # Get the number of models in the file:
    entities = gmsh.model.getEntities()

    # Get the number of physical groups in the file:
    physical_groups = gmsh.model.getPhysicalGroups()

    # Dictionary to store the count of elements by dimension
    nodes_by_dim = {}
    elements_by_dim = {}  

    print("\nComplete Info: \n---------------------------------")    

    for e in entities:
        # Dimension and tag of the entity:
        dim = e[0]
        tag = e[1]

        # Get the mesh nodes for the entity (dim, tag):
        nodeTags, nodeCoords, nodeParams = gmsh.model.mesh.getNodes(dim, tag)
        if dim not in nodes_by_dim:
            nodes_by_dim[dim] = 0
        for tags in nodeTags:
            nodes_by_dim[dim] += 1

        # Get the mesh elements for the entity (dim, tag):
        elemTypes, elemTags, elemNodeTags = gmsh.model.mesh.getElements(dim, tag)
        if dim not in elements_by_dim:
            elements_by_dim[dim] = 0
        for tags in elemTags:
            elements_by_dim[dim] += len(tags)  

        # * Type and name of the entity:
        info_mode = gmsh.model.getType(dim, tag)
        name = gmsh.model.getEntityName(dim, tag)
        if len(name): name += ' '
        print("Entity " + name + str(e) + " of type " + info_mode)

        # * Number of mesh nodes and elements:
        numElem = sum(len(i) for i in elemTags)
        print(" - Mesh has " + str(len(nodeTags)) + " nodes and " + str(numElem) +
            " elements")

        # * Upward and downward adjacencies:
        up, down = gmsh.model.getAdjacencies(dim, tag)
        if len(up):
            print(" - Upward adjacencies: " + str(up))
        if len(down):
            print(" - Downward adjacencies: " + str(down))

        # * Does the entity belong to physical groups?
        physicalTags = gmsh.model.getPhysicalGroupsForEntity(dim, tag)
        if len(physicalTags):
            s = ''
            for p in physicalTags:
                n = gmsh.model.getPhysicalName(dim, p)
                if n: n += ' '
                s += n + '(' + str(dim) + ', ' + str(p) + ') '
            print(" - Physical groups: " + s)

        # * Is the entity a partition entity? If so, what is its parent entity?
        partitions = gmsh.model.getPartitions(dim, tag)
        if len(partitions):
            print(" - Partition tags: " + str(partitions) + " - parent entity " +
                str(gmsh.model.getParent(dim, tag)))

        # * List all types of elements making up the mesh of the entity:
        for t in elemTypes:
            name, dim, order, numv, parv, _ = gmsh.model.mesh.getElementProperties(t)
            print(" - Element type: " + name + ", order " + str(order) + " (" +
                str(numv) + " nodes in param coord: " + str(parv) + ")") 
            
    # Display node count by dimension
    for dim, count in nodes_by_dim.items():
        print(f"Resume     : {count} nodes in dimension {dim}")

    # Display element count by dimension
    for dim, count in elements_by_dim.items():
        print(f"Resume     : {count} elements in dimension {dim}")

# %% [markdown]
# # `get_data()`

# %%
def get_data(FINITE_ELEMENT, BOUNDARY, MATERIAL, model, info_mode=False):
    element, order = FINITE_ELEMENT
    file_path = f"pre_processing/mesh/{model}_domain_{element}{order}.msh"
    
    gmsh.initialize()
    gmsh.open(file_path)
    basic_info()
    if info_mode:
        complete_info()

    # Structure Data
    mesh_data = {}
    mesh_data['cell'] = get_cell_data(MATERIAL)
    mesh_data['nodes'] = get_nodes_data(BOUNDARY)
    mesh_data['edges'] = get_edge_data()

    gmsh.finalize()    
    return mesh_data

# %% [markdown]
# # `get_tetrahedra_data()`

# %%
def get_tetrahedra_data(FINITE_ELEMENT, BOUNDARY, MATERIAL, model):    
    gmsh.initialize()
    element, order = FINITE_ELEMENT
    gmsh.open(f"pre_processing/mesh/{model}_domain_{element}{order}.msh")
    basic_info()

    # Structure Data
    mesh_data = {}
    mesh_data['cell'] = get_new_cell_data(MATERIAL, problem_dim=3)
    mesh_data['nodes'] = get_nodes_data(BOUNDARY, problem_dim=2)
    # mesh_data['edges'] = get_edge_data()

    gmsh.finalize()    
    return mesh_data

# %% [markdown]
# Conversão do arquivo Jupyter Notebook para um script Python: ``python -m nbconvert --to script name.ipynb``
# 
# Belo Horizonte, Brazil.  
# Adilton Junio Ladeira Pereira - adt@ufmg.br  
# &copy; All rights reserved.
# 
# version 1.0. November, 2024.


