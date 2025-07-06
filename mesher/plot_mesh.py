import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from matplotlib.patches import Patch

def structured_data(mesh_data):
    # Extraindo as coordenadas globais dos nós (x, y)
    try:
        xg = [node['xg'][0] for node in mesh_data['nodes'].values()]
        yg = [node['xg'][1] for node in mesh_data['nodes'].values()]
    except KeyError as e:
        raise ValueError(f"Erro ao acessar as coordenadas globais: {e}")

    # Extraindo a matriz de conectividade
    try:
        conn = [cell['conn'] for cell in mesh_data['cell'].values()]
        # Ajusta índice para 0-based
        conn_py = [[node - 1 for node in nodes[:3]] for nodes in conn]  
    except KeyError as e:
        raise ValueError(f"Erro ao acessar a conectividade: {e}")
    
    return xg, yg, conn_py


def plot_triangular_mesh(INFO_GRAPH, mesh_data, nodes_index_based=0):
    # Dados do gráfico
    show_cell = INFO_GRAPH['cell']
    show_nodes = INFO_GRAPH['nodes']
    show_edges = INFO_GRAPH['edges']
    show_edges_numb = INFO_GRAPH['edges_numb']
    filepath = INFO_GRAPH['filepath']

    # Estruturando os dados da malha
    nodes_data = mesh_data['nodes']
    
    # Extraindo as coordenadas globais dos nós (x, y) e a matriz de conectividade
    xg, yg, conn_py = structured_data(mesh_data)

    # Plotando a malha de elementos finitos
    plt.figure(figsize=(8, 6))
    plt.triplot(xg, yg, conn_py, color='gray')  
    
    # Plotando as arestas dos elementos
    if show_edges or show_edges_numb:
        for key, edge in mesh_data['edges'].items():
            # Coordenadas dos nós inicial e final
            x0 = nodes_data[edge['conn'][0]]['xg'][0]
            y0 = nodes_data[edge['conn'][0]]['xg'][1]
            x1 = nodes_data[edge['conn'][1]]['xg'][0]
            y1 = nodes_data[edge['conn'][1]]['xg'][1]
            
            # Ponto médio da aresta
            x_mid, y_mid = (x0 + x1) / 2, (y0 + y1) / 2
            
            # Vetor da seta (a partir do ponto médio)
            dx, dy = (x1 - x0) * 0.2, (y1 - y0) * 0.2  
                
            # Adicionando uma seta no meio da aresta
            if show_edges:
                plt.arrow(x_mid, y_mid, dx, dy, head_width=0.015, head_length=0.05,
                                fc='blue', ec='blue', length_includes_head=True)

            # Adicionando os números das arestas
            if show_edges_numb:
                plt.scatter(x_mid, y_mid, marker='s', color='white', 
                                edgecolor='black', s=120, zorder=1)                
                plt.text(x_mid, y_mid, key, color='blue', fontsize=6, 
                                ha='center', va='center')

    # Adicionando nós
    if show_nodes: 
        for key, node in nodes_data.items():
            x, y = node['xg'][0], node['xg'][1]
            plt.scatter(x, y, color='white', edgecolor='black', s=180)
            if nodes_index_based == 1:
                plt.text(x, y, str(key), color='red', fontsize=6, ha='center', va='center')
            elif nodes_index_based == 0:
                plt.text(x, y, str(key - 1), color='red', fontsize=6, ha='center', va='center')
            else:
                raise ValueError("nodes_index_based deve ser 0 ou 1.")  
    else:
        plt.scatter(xg, yg, color='black', s=1, zorder=3)

    # Adicionando elementos
    if show_cell:
        for key, cell in mesh_data['cell'].items():
            xc = np.mean([nodes_data[node]['xg'][0] for node in cell['conn']])
            yc = np.mean([nodes_data[node]['xg'][1] for node in cell['conn']])
            plt.text(xc, yc, str(key), fontweight='bold',
                        color='black', fontsize=6, ha='center', va='center')
                
    # Ajustando rótulos e layout
    plt.xlabel(r'$x$')
    plt.ylabel(r'$y$')
    plt.axis('equal')
    plt.tight_layout()
    
    # Salvando o arquivo no formato SVG
    plt.savefig(filepath, format="svg")
    plt.show()
    print(f"Arquivo gráfico da malha salvo em: {filepath}")


def get_faces(elem_vertices):
    """
    Retorna as 4 faces de um tetraedro de acordo com a convenção
    deduzida da estrutura de dados.
    elem_vertices: uma lista ou array com os 4 índices dos vértices do elemento, ex: [v0, v1, v2, v3]
    """
    v0, v1, v2, v3 = elem_vertices
    
    # Convenção deduzida:
    # Face 0 é oposta ao v3
    # Face 1 é oposta ao v2
    # Face 2 é oposta ao v0
    # Face 3 é oposta ao v1
    faces = [
        (v0, v1, v2),  # Face 0
        (v0, v1, v3),  # Face 1
        (v1, v2, v3),  # Face 2
        (v0, v2, v3)   # Face 3
    ]
    return faces


def get_edges(EToV):
    """
    Extrai todas as arestas únicas da malha para desenhar o contorno.
    """
    all_edges = set()
    for elem in EToV:
        # Usamos a mesma convenção de faces para extrair as arestas
        faces = get_faces(elem)
        for face in faces:
            # Adiciona arestas da face em ordem canônica (menor, maior) para evitar duplicatas
            all_edges.add(tuple(sorted((face[0], face[1]))))
            all_edges.add(tuple(sorted((face[1], face[2]))))
            all_edges.add(tuple(sorted((face[2], face[0]))))
    return list(all_edges)


def plot_cubeK6_mesh(VX, VY, VZ, EToV):
    # Arestas do cubo (pares de índices dos nós)
    cube_edges = get_edges(EToV)

    fig = plt.figure(figsize=(10, 8))
    colors = plt.cm.tab10(np.linspace(0, 1, len(EToV)))
    for k, elem in enumerate(EToV):
        ax = fig.add_subplot(2, 3, k + 1, projection='3d')
        # ax.view_init(elev=0, azim=30)

        # 1. Desenhar arestas do cubo (tracejadas)
        for i, j in cube_edges:
            x = [VX[i], VX[j]]
            y = [VY[i], VY[j]]
            z = [VZ[i], VZ[j]]
            ax.plot(x, y, z, 'k--', linewidth=0.6, zorder=0)

        # 2. Desenhar tetraedro colorido
        faces = get_faces(elem)
        for face in faces:
            verts = [[VX[i], VY[i], VZ[i]] for i in face]
            poly = Poly3DCollection([verts], facecolor='b', edgecolor='k', alpha=0.7)
            ax.add_collection3d(poly)

        # 3. Adicionar vértices e rótulos dos nós
        ax.scatter(VX[elem], VY[elem], VZ[elem], color='black', s=20)
        for idx in elem:
            ax.text(VX[idx]*1.1, VY[idx]*1.1, VZ[idx]*1.1, str(idx), color='red',
                    ha='center', va='center', fontsize=10, weight='bold', zorder=1)

        # 4. Rótulo do ID do elemento
        centroid = np.mean([[VX[i], VY[i], VZ[i]] for i in elem], axis=0)
        ax.text(*centroid, f't{k}', color='white', ha='center', va='center',
                fontsize=14, weight='bold',
                bbox=dict(facecolor='black', alpha=0.5, boxstyle='circle'))

        # 5. Ajuste dos eixos
        ax.set_xlim(-0.1, 1.1); ax.set_ylim(-0.1, 1.1); ax.set_zlim(-0.1, 1.1)
        ax.set_box_aspect([1, 1, 1])
        ax.axis('off')

    fig.suptitle("Unit cube", fontsize=12)
    plt.tight_layout()


def plot_local_cubeK6_mesh(VX, VY, VZ, EToV):
    """
    Plota a malha tetraédrica, ordenando as faces de cada elemento
    de trás para a frente para evitar sobreposição incorreta de cores (Z-fighting).
    """
    mesh_edges = get_edges(EToV)
    face_colors = ['red', 'green', 'blue', 'gold']
    face_labels = ['Face 0 (oposta a v3)', 'Face 1 (oposta a v2)', 'Face 2 (oposta a v0)', 'Face 3 (oposta a v1)']

    fig = plt.figure(figsize=(12, 8))
    
    for k, elem in enumerate(EToV):
        ax = fig.add_subplot(2, 3, k + 1, projection='3d')
        
        # Define um ângulo de visão fixo para garantir consistência
        # ax.view_init(elev=25, azim=45)

        # 1. Desenhar arestas de contorno
        for i, j in mesh_edges:
            x, y, z = ([VX[i], VX[j]], [VY[i], VY[j]], [VZ[i], VZ[j]])
            ax.plot(x, y, z, 'k--', linewidth=0.8, alpha=0.5, zorder=0)

        # 2. Lógica para ordenar e desenhar as faces (O "Algoritmo do Pintor")
        faces = get_faces(elem)
        
        # Obter o vetor de visão da câmera a partir dos ângulos do subplot
        elev = np.deg2rad(ax.elev)
        azim = np.deg2rad(ax.azim)
        view_vector = np.array([np.cos(elev) * np.cos(azim),
                                np.cos(elev) * np.sin(azim),
                                np.sin(elev)])
        
        # Calcular a profundidade de cada face e armazená-la
        face_info = []
        for face_idx, face in enumerate(faces):
            # Calcular o centroide da face
            centroid = np.mean([[VX[i], VY[i], VZ[i]] for i in face], axis=0)
            # Calcular a 'profundidade' como o produto escalar com o vetor de visão
            # Um valor menor significa que está "mais longe" da câmera
            depth = np.dot(centroid, view_vector)
            face_info.append({'face': face, 'idx': face_idx, 'depth': depth})
        
        # Ordenar as faces pela profundidade, da mais distante para a mais próxima
        sorted_faces = sorted(face_info, key=lambda x: x['depth'])

        # Desenhar as faces na ordem correta (de trás para a frente)
        for info in sorted_faces:
            face = info['face']
            face_idx = info['idx']
            verts = [[VX[i], VY[i], VZ[i]] for i in face]
            poly = Poly3DCollection([verts], facecolor=face_colors[face_idx], edgecolor='k', alpha=0.8)
            ax.add_collection3d(poly)

        # 3. Adicionar vértices e rótulos dos nós
        ax.scatter(VX[elem], VY[elem], VZ[elem], color='black', s=50)
        for idx in elem:
            ax.text(VX[idx]*1.1, VY[idx]*1.1, VZ[idx]*1.1, str(idx), color='black',
                    ha='center', va='center', fontsize=10, weight='bold')

        # 4. Rótulo do ID do elemento
        centroid = np.mean([[VX[i], VY[i], VZ[i]] for i in elem], axis=0)
        ax.text(*centroid, f't{k}', color='white', ha='center', va='center',
                fontsize=14, weight='bold',
                bbox=dict(facecolor='black', alpha=0.5, boxstyle='circle'))

        # 5. Ajuste dos eixos
        ax.set_xlim(-0.1, 1.1); ax.set_ylim(-0.1, 1.1); ax.set_zlim(-0.1, 1.1)
        ax.set_box_aspect([1, 1, 1])
        ax.axis('off')

    legend_elements = [Patch(facecolor=color, edgecolor='k', label=label)
                       for color, label in zip(face_colors, face_labels)]
    fig.axes[0].legend(handles=legend_elements, loc='upper left', bbox_to_anchor=(-0.1, 1.1))

    fig.suptitle("Unit Cubic with local faces", fontsize=12)
    plt.tight_layout(rect=[0, 0, 1, 0.96])


def extract_all_mesh_edges(EToV):
    """
    Extrai todas as arestas únicas de uma malha tetraédrica.

    Parâmetros:
        EToV (ndarray): Matriz de conectividade (Elementos x Vértices).

    Retorna:
        list: Uma lista de tuplas, onde cada tupla (i, j) representa uma aresta única.
    """
    # Um set é usado para armazenar as arestas e garantir a unicidade automaticamente.
    edges_set = set()
    
    # Itera sobre cada elemento (tetraedro)
    for element_vertices in EToV:
        # Um tetraedro com vértices v0, v1, v2, v3 tem 6 arestas
        v0, v1, v2, v3 = element_vertices
        element_edges = [
            (v0, v1), (v0, v2), (v0, v3),
            (v1, v2), (v1, v3),
            (v2, v3)
        ]
        
        # Adiciona as arestas ao set em um formato canônico (ordenado)
        # para evitar duplicatas como (a, b) e (b, a).
        for i, j in element_edges:
            edges_set.add(tuple(sorted((i, j))))
            
    return list(edges_set)


def plot_wireframe(ax, VX, VY, VZ, edges):
    """
    Desenha um wireframe (conjunto de arestas) em um eixo 3D.
    Esta função agora desenha o wireframe completo da malha.
    """
    for i, j in edges:
        x = [VX[i], VX[j]]
        y = [VY[i], VY[j]]
        z = [VZ[i], VZ[j]]
        ax.plot(x, y, z, 'k--', linewidth=0.6, alpha=0.5, zorder=0)


def get_domain_cube_geometry(domain_limits=((-0.1, 1.1), (-0.1, 1.1), (-0.1, 1.1))):
    """
    Cria a geometria de um cubo de domínio para referência.

    Parâmetros:
        domain_limits (tuple): Tupla contendo os limites (min, max) para x, y, e z.

    Retorna:
        tuple: (vértices do domínio, arestas do domínio)
    """
    xmin, xmax = domain_limits[0]
    ymin, ymax = domain_limits[1]
    zmin, zmax = domain_limits[2]

    # 8 vértices do cubo de domínio
    verts = np.array([
        [xmin, ymin, zmin], [xmax, ymin, zmin], [xmin, ymax, zmin], [xmax, ymax, zmin],
        [xmin, ymin, zmax], [xmax, ymin, zmax], [xmin, ymax, zmax], [xmax, ymax, zmax]
    ])
    # 12 arestas que conectam os 8 vértices
    edges = [
        (0, 1), (0, 2), (1, 3), (2, 3), (4, 5), (4, 6),
        (5, 7), (6, 7), (0, 4), (1, 5), (2, 6), (3, 7)
    ]
    return verts, edges


def plot_element_faces(ax, element_vertices):
    """Desenha as faces de um elemento tetraédrico."""
    # Um tetraedro tem 4 faces triangulares
    faces_indices = [[0, 1, 2], [0, 1, 3], [0, 2, 3], [1, 2, 3]]
    
    # Coleta as coordenadas dos vértices para cada face
    faces_verts = [[element_vertices[i] for i in face] for face in faces_indices]

    poly = Poly3DCollection(faces_verts, facecolor='cyan', edgecolor='k', alpha=0.2)
    ax.add_collection3d(poly)


def plot_dg_nodes(ax, node_coords, node_ids):
    """Plota os nós de DG e seus IDs globais."""
    ax.scatter(node_coords[:, 0], node_coords[:, 1], node_coords[:, 2], color='red', s=18, zorder=10)
    for i in range(len(node_ids)):
        ax.text(node_coords[i, 0], node_coords[i, 1], node_coords[i, 2] + 0.08,
                f'{node_ids[i]}', color='black', ha='center', va='bottom', fontsize=8, weight='bold')


def plot_element_label(ax, element_vertices, element_id):
    """Adiciona um rótulo com o ID do elemento no seu centroide."""
    # 1. Encontrar o vértice com a maior coordenada Z
    top_vertex_index = np.argmax(element_vertices[:, 2])
    top_vertex_coords = element_vertices[top_vertex_index]

    # 2. Posicionar o texto um pouco acima deste vértice
    text_position = top_vertex_coords + np.array([0, 0, 1.0]) # Offset em Z

    ax.text(*text_position, f't{element_id}', color='white', ha='center', va='center',
            fontsize=8, weight='bold', bbox=dict(facecolor='navy', alpha=0.6, boxstyle='circle'))


def setup_plot_axes(ax, limits):
    """Configura os limites, proporção e aparência dos eixos."""
    ax.set_xlim(limits[0])
    ax.set_ylim(limits[1])
    ax.set_zlim(limits[2])

    # Correção: Garante que o aspect ratio seja um vetor 1D de shape (3,)
    aspect_ratio = np.array([
        np.diff(limits[0]),
        np.diff(limits[1]),
        np.diff(limits[2])
    ]).flatten()  # O método .flatten() converte o array (3, 1) para (3,)

    ax.set_box_aspect(aspect_ratio)
    ax.axis('off')


def find_interface_elements_on_face(dg_solver, group_id, face_index, group_cube_dim):
    """
    Encontra os elementos de um grupo que estão em uma face de interface específica.

    Esta função encapsula a lógica de filtragem para identificar quais
    elementos de um 'group_id' tocam uma face específica do cubo que o define.

    Retorna:
        list: Uma lista de índices de elementos que satisfazem os critérios.
    """
    if group_cube_dim is None or group_cube_dim <= 0:
        raise ValueError("A dimensão do cubo do grupo físico (group_cube_dim) deve ser um número positivo.")

    # Mapeia o índice da face (0-5) para um eixo e direção
    face_definitions = {
        0: {'axis': 0, 'sign': +1, 'name': '+X'}, 1: {'axis': 0, 'sign': -1, 'name': '-X'},
        2: {'axis': 1, 'sign': +1, 'name': '+Y'}, 3: {'axis': 1, 'sign': -1, 'name': '-Y'},
        4: {'axis': 2, 'sign': +1, 'name': '+Z'}, 5: {'axis': 2, 'sign': -1, 'name': '-Z'}
    }
    if face_index not in face_definitions:
        raise ValueError(f"Índice de face (No_FACE) inválido: {face_index}. Deve ser entre 0 e 5.")
    
    face_def = face_definitions[face_index]
    plane_coord = face_def['sign'] * (group_cube_dim / 2.0)
    plane_axis = face_def['axis']

    elements_on_face = []
    target_elements = dg_solver.mesh.get_elements_by_group(group_id)

    # Itera sobre os elementos do grupo para encontrar os que estão na interface correta
    for k in target_elements:
        for f in range(dg_solver.n_faces):
            k_neighbor = dg_solver.EToE[k, f]
            # É uma face de interface?
            if k != k_neighbor and dg_solver.mesh.EToG[k_neighbor] != group_id:
                # Calcular o centroide desta face para verificar sua posição
                face_nodes_local = dg_solver.fmask[:, f]
                coords_k = np.array([dg_solver.x[:, k], dg_solver.y[:, k], dg_solver.z[:, k]]).T
                face_centroid = np.mean(coords_k[face_nodes_local], axis=0)
                
                # A face está no plano correto?
                if np.isclose(face_centroid[plane_axis], plane_coord, atol=1e-9):
                    elements_on_face.append(k)
                    break  # Elemento encontrado, otimiza o laço

    return sorted(list(set(elements_on_face)))


def get_group_wireframe_edges(dg_solver, group_id):
    """
    Prepara as arestas do wireframe para um grupo físico específico.
    """
    # Obter os elementos que pertencem ao grupo de interesse
    group_element_indices = dg_solver.mesh.get_elements_by_group(group_id)
    if len(group_element_indices) == 0:
        return []
    
    # Criar uma matriz de conectividade contendo apenas esses elementos
    group_EToV = dg_solver.mesh.EToV[group_element_indices]
    
    # Extrair as arestas apenas deste subconjunto de elementos
    return extract_all_mesh_edges(group_EToV)


def plot_buildMaps_cubeK6(dg_solver):
    """
    Plota a malha e os nós de DG, usando o wireframe completo da malha como referência.
    Esta versão implementa o comportamento visual solicitado.
    """
    # Extrai os dados básicos do solver
    EToV = dg_solver.mesh.EToV
    VX, VY, VZ = dg_solver.mesh.vx, dg_solver.mesh.vy, dg_solver.mesh.vz
    K = dg_solver.mesh.number_of_elements()

    # 1. Extrai todas as arestas únicas da malha para o wireframe
    all_mesh_edges = extract_all_mesh_edges(EToV)
    
    # Define os limites do plot com base na extensão total dos vértices
    domain_limits = (
        (VX.min() - 0.1, VX.max() + 0.1),
        (VY.min() - 0.1, VY.max() + 0.1),
        (VZ.min() - 0.1, VZ.max() + 0.1)
    )
    
    # Cria a figura principal
    fig = plt.figure(figsize=(18, 10)) # Aumentei o tamanho para melhor visualização
    fig.suptitle("Malha com Índices Globais dos Nós (Comportamento Corrigido)", fontsize=16)

    # Itera sobre cada elemento para criar um subplot
    for k in range(K):
        ax = fig.add_subplot(2, 3, k + 1, projection='3d')

        # Coleta de dados específicos do elemento 'k'
        vert_indices = EToV[k]
        element_vertices = np.array([VX[vert_indices], VY[vert_indices], VZ[vert_indices]]).T
        node_coords_k = np.array([dg_solver.x[:, k], dg_solver.y[:, k], dg_solver.z[:, k]]).T
        node_ids_k = dg_solver.node_ids[:, k]

        # --- Orquestração da plotagem ---
        # 2. Desenha o wireframe completo da malha como fundo
        plot_wireframe(ax, VX, VY, VZ, all_mesh_edges)
        
        # 3. Desenha o elemento atual em destaque
        plot_element_faces(ax, element_vertices)
        
        # 4. Desenha os nós de DG e seus rótulos
        plot_dg_nodes(ax, node_coords_k, node_ids_k)
        
        # 5. Adiciona o rótulo do elemento
        plot_element_label(ax, element_vertices, k)
        
        # 6. Configura os eixos
        setup_plot_axes(ax, domain_limits)

    plt.tight_layout(rect=[0, 0, 1, 0.96])


def plot_buildInterfaceMaps_cubeK96(dg_solver, GROUP_ID=1, No_FACE=0, tfz_dim=None):
    """
    Plota os elementos de uma interface específica usando funções auxiliares refatoradas.
    """
    # 1. Encontrar os elementos de interesse usando a nova função de filtragem
    elements_to_plot = find_interface_elements_on_face(dg_solver, GROUP_ID, No_FACE, tfz_dim)
    
    if not elements_to_plot:
        print(f"Nenhum elemento do grupo {GROUP_ID} encontrado na face {No_FACE}.")
        return

    print(f"\n Encontrados {len(elements_to_plot)} elementos na interface. Plotando...")

    # 2. Preparar o wireframe de fundo para o grupo de interesse
    background_edges = get_group_wireframe_edges(dg_solver, GROUP_ID)

    # 3. Preparar dados gerais para plotagem
    VX, VY, VZ = dg_solver.mesh.vx, dg_solver.mesh.vy, dg_solver.mesh.vz
    domain_limits = ((VX.min()-0.1, VX.max()+0.1), (VY.min()-0.1, VY.max()+0.1), (VZ.min()-0.1, VZ.max()+0.1))
    
    # 4. Configurar a figura
    face_name = {0:'+X', 1:'-X', 2:'+Y', 3:'-Y', 4:'+Z', 5:'-Z'}.get(No_FACE, '')
    fig, axes = plt.subplots(2, 2, subplot_kw={'projection': '3d'}, figsize=(12, 10))
    axes = axes.flatten()
    fig.suptitle(f"Elementos do Grupo {GROUP_ID} na Interface da Face '{face_name}'", fontsize=12)

    # 5. Loop de plotagem (orquestração)
    for i, k in enumerate(elements_to_plot):
        if i >= 4: break
        ax = axes[i]
        
        # Coleta de dados do elemento 'k'

        # 1. Obtém os índices dos 4 vértices que formam o elemento 'k'
        vert_indices = dg_solver.mesh.EToV[k]

        # 2. Usa os índices para buscar as coordenadas (x,y,z) desses vértices
        element_vertices = np.array([VX[vert_indices], VY[vert_indices], VZ[vert_indices]]).T

        # 3. Coleta as coordenadas e IDs de todos os nós de DG dentro do elemento 'k'
        node_coords_k = np.array([dg_solver.x[:, k], dg_solver.y[:, k], dg_solver.z[:, k]]).T
        node_ids_k = dg_solver.node_ids[:, k]

        # Chamada às funções de plotagem modulares
        plot_wireframe(ax, VX, VY, VZ, background_edges)
        plot_element_faces(ax, element_vertices)
        plot_dg_nodes(ax, node_coords_k, node_ids_k)
        plot_element_label(ax, element_vertices, k)
        setup_plot_axes(ax, domain_limits)
    
    # Limpa eixos não utilizados
    for j in range(len(elements_to_plot), 4):
        axes[j].axis('off')
        
    plt.tight_layout(rect=[0, 0, 1, 0.95])

