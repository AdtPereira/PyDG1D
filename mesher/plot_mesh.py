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


def plot_dg_nodes(ax, node_coords, node_ids, highlight_nodes=None):
    """
    Plota os nós de DG e seus IDs globais, com a opção de destacar um subconjunto.

    Parâmetros:
        ax: Eixo do Matplotlib.
        node_coords (ndarray): Coordenadas (x,y,z) dos nós do elemento.
        node_ids (ndarray): IDs globais dos nós.
        highlight_nodes (set, optional): Um set de IDs de nós a serem destacados.
    """
    if highlight_nodes is None:
        highlight_nodes = set()

    # Plota os nós de DG e seus IDs globais
    for i in range(len(node_ids)):
        node_id = node_ids[i]
        coords = node_coords[i]
        
        if node_id in highlight_nodes:
            # Nó destacado (de fronteira)
            ax.scatter(coords[0], coords[1], coords[2], color='magenta', s=60, zorder=12, edgecolor='black')
        else:
            # Nó padrão (interno ao elemento)
            ax.scatter(coords[0], coords[1], coords[2], color='red', s=20, zorder=10)
            
        ax.text(coords[0], coords[1], coords[2] + 0.08,
                f'{node_id}', color='black', ha='center', va='bottom', fontsize=8, weight='bold')


def plot_element_label(ax, element_vertices, element_id):
    """Adiciona um rótulo com o ID do elemento no seu centroide."""
    # 1. Encontrar o vértice com a maior coordenada Z
    top_vertex_index = np.argmax(element_vertices[:, 2])
    top_vertex_coords = element_vertices[top_vertex_index]

    # 2. Posicionar o texto um pouco acima deste vértice
    text_position = top_vertex_coords + np.array([-0.8, 0, 0.0]) # Offset em Z
    
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


def get_group_wireframe_edges(dg, group_id):
    """
    Prepara as arestas do wireframe para um grupo físico específico.
    """
    # Obter os elementos que pertencem ao grupo de interesse
    group_element_indices = dg.mesh.get_elements_by_group(group_id)
    if len(group_element_indices) == 0:
        return []
    
    # Criar uma matriz de conectividade contendo apenas esses elementos
    group_EToV = dg.mesh.EToV[group_element_indices]
    
    # Extrair as arestas apenas deste subconjunto de elementos
    return extract_all_mesh_edges(group_EToV)


def plot_buildMaps_cubeK6(dg):
    """
    Plota a malha e os nós de DG, usando o wireframe completo da malha como referência.
    Esta versão implementa o comportamento visual solicitado.
    """
    # Extrai os dados básicos do solver
    EToV = dg.mesh.EToV
    VX, VY, VZ = dg.mesh.vx, dg.mesh.vy, dg.mesh.vz
    K = dg.mesh.number_of_elements()

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
        node_coords_k = np.array([dg.x[:, k], dg.y[:, k], dg.z[:, k]]).T
        node_ids_k = dg.node_ids[:, k]

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


def plot_buildInterfaceMaps_cubeK96(dg, GROUP_ID=1, No_FACE=0, tfz_dim=None):
    """
    Plota os elementos de uma interface específica usando funções auxiliares refatoradas.
    """
    # 1. Encontrar os elementos de interesse usando a nova função de filtragem
    elements_to_plot = find_interface_elements_on_face(dg, GROUP_ID, No_FACE, tfz_dim)
    
    if not elements_to_plot:
        print(f"Nenhum elemento do grupo {GROUP_ID} encontrado na face {No_FACE}.")
        return

    print(f"\n Encontrados {len(elements_to_plot)} elementos na interface. Plotando...")

    # 2. Preparar o wireframe de fundo para o grupo de interesse
    background_edges = get_group_wireframe_edges(dg, GROUP_ID)

    # 3. Preparar dados gerais para plotagem
    VX, VY, VZ = dg.mesh.vx, dg.mesh.vy, dg.mesh.vz
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
        vert_indices = dg.mesh.EToV[k]

        # 2. Usa os índices para buscar as coordenadas (x,y,z) desses vértices
        element_vertices = np.array([VX[vert_indices], VY[vert_indices], VZ[vert_indices]]).T

        # 3. Coleta as coordenadas e IDs de todos os nós de DG dentro do elemento 'k'
        node_coords_k = np.array([dg.x[:, k], dg.y[:, k], dg.z[:, k]]).T
        node_ids_k = dg.node_ids[:, k]

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


def validate_plot_arguments(plot_type, layout):
    """Valida os argumentos de tipo de plotagem e layout."""
    valid_plot_types = ['all', 'boundary', 'internal']
    if plot_type not in valid_plot_types:
        raise ValueError(f"Tipo de plotagem inválido '{plot_type}'. Use um de {valid_plot_types}.")
    
    valid_layouts = ['separated', 'combined']
    if layout not in valid_layouts:
        raise ValueError(f"Layout inválido '{layout}'. Use 'separated' ou 'combined'.")


def classify_elements(dg):
    """Classifica os elementos da malha em internos e de fronteira."""
    all_elements_set = set(range(dg.mesh.number_of_elements()))
    
    boundary_elements_set = set()
    boundary_node_ids = set()
    
    # Verifica se existem nós de fronteira definidos
    if dg.mapB.size > 0:
        # Calcula os índices dos elementos que contêm os nós de fronteira
        elem_indices = dg.mapB // (dg.n_fp * dg.n_faces)
        boundary_elements_set = set(np.unique(elem_indices))
        boundary_node_ids = set(dg.vmapB)
    
    internal_elements_set = all_elements_set - boundary_elements_set
    
    return boundary_elements_set, internal_elements_set, boundary_node_ids


def select_elements_to_plot(plot_type, boundary_set, internal_set):
    """Seleciona os elementos a serem plotados e define o título e a legenda."""
    if plot_type == 'all':
        elements = sorted(list(boundary_set | internal_set))
        title = "Visualização de Todos os Elementos"
        patches = [
            Patch(facecolor='cyan', label='Elemento de Fronteira'),
            Patch(facecolor='lightgreen', label='Elemento Interno')
        ]
    elif plot_type == 'boundary':
        elements = sorted(list(boundary_set))
        title = "Visualização dos Elementos de Fronteira"
        patches = [Patch(facecolor='cyan', label='Elemento de Fronteira')]
    else:  # 'internal'
        elements = sorted(list(internal_set))
        title = "Visualização dos Elementos Internos"
        patches = [Patch(facecolor='lightgreen', label='Elemento Interno')]
        
    return elements, title, patches


def setup_figure(layout, num_elements, title):
    """Cria e configura a figura e os eixos do Matplotlib."""
    if layout == 'separated':
        ncols = min(num_elements, 4)
        nrows = (num_elements + ncols - 1) // ncols
        fig, axes = plt.subplots(nrows, ncols, subplot_kw={'projection': '3d'}, figsize=(18, 4.5 * nrows))
        axes = axes.flatten() if num_elements > 1 else [axes]
    else:  # 'combined'
        fig = plt.figure(figsize=(12, 10))
        axes = [fig.add_subplot(111, projection='3d')]
        
    fig.suptitle(title, fontsize=16)
    return fig, axes


def draw_single_element(ax, dg, k, plot_config):
    """Desenha um único elemento (faces, nós, rótulos) em um eixo específico."""
    # Extrai dados do elemento 'k'
    vert_indices = dg.mesh.EToV[k]
    element_vertices = np.array([dg.mesh.vx[vert_indices], dg.mesh.vy[vert_indices], dg.mesh.vz[vert_indices]]).T
    node_coords_k = np.array([dg.x[:, k], dg.y[:, k], dg.z[:, k]]).T
    node_ids_k = dg.node_ids[:, k]

    # Determina a cor e quais nós destacar
    if k in plot_config['boundary_elements_set']:
        face_color, nodes_to_highlight = 'cyan', plot_config['boundary_node_ids']
    else:
        face_color, nodes_to_highlight = 'lightgreen', None

    # Desenha o wireframe de fundo
    plot_wireframe(ax, dg.mesh.vx, dg.mesh.vy, dg.mesh.vz, plot_config['all_mesh_edges'])

    # Desenha as faces do elemento
    faces_indices = [[0, 1, 2], [0, 1, 3], [0, 2, 3], [1, 2, 3]]
    faces_verts = [[element_vertices[i] for i in face] for face in faces_indices]
    poly = Poly3DCollection(faces_verts, facecolor=face_color, edgecolor='k', alpha=plot_config['alpha'])
    ax.add_collection3d(poly)

    # Desenha os nós de DG se aplicável
    if plot_config['should_plot_nodes']:
        plot_dg_nodes(ax, node_coords_k, node_ids_k, highlight_nodes=nodes_to_highlight)

    # Adiciona rótulos e configura os eixos
    plot_element_label(ax, element_vertices, k)
    setup_plot_axes(ax, plot_config['domain_limits'])
    if plot_config['layout'] == 'separated':
        ax.set_title(f"Elemento {k}", fontsize=10)


def plot_element_faces(ax, element_vertices, face_color='cyan', alpha=0.2):
    """Desenha as faces de um elemento tetraédrico."""
    faces_indices = [[0, 1, 2], [0, 1, 3], [0, 2, 3], [1, 2, 3]]
    faces_verts = [[element_vertices[i] for i in face] for face in faces_indices]
    poly = Poly3DCollection(faces_verts, facecolor=face_color, edgecolor='k', alpha=alpha)
    ax.add_collection3d(poly)


def plot_elements(dg, plot_type='all', layout='separated', ELEMENT_ID=None):
    """
    Plota os elementos da malha com opções de filtragem e layout.

    Parâmetros:
        dg: Objeto da classe Maxwell3D.
        plot_type (str): 'all', 'boundary' ou 'internal'. Ignorado se ELEMENT_ID for fornecido.
        layout (str): 'separated' ou 'combined'.
        ELEMENT_ID (int ou list, opcional): Um ID de elemento específico ou uma lista de IDs
                                             para plotar. Se fornecido, sobrepõe 'plot_type'.
    """
    # 1. Validação dos argumentos de entrada
    validate_plot_arguments(plot_type, layout)

    # 2. Classificação de todos os elementos (necessário para cores e legendas)
    boundary_elements, internal_elements, boundary_nodes = classify_elements(dg)

    # 3. Seleção de elementos com base em ELEMENT_ID ou plot_type
    if ELEMENT_ID is not None:
        # Padroniza ELEMENT_ID para ser sempre uma lista
        if isinstance(ELEMENT_ID, int):
            elements_to_plot = [ELEMENT_ID]
        elif isinstance(ELEMENT_ID, list):
            elements_to_plot = ELEMENT_ID
        else:
            raise TypeError("ELEMENT_ID deve ser um inteiro ou uma lista de inteiros.")

        # Valida se os IDs fornecidos são válidos
        max_id = dg.mesh.number_of_elements() - 1
        if any(not (0 <= eid <= max_id) for eid in elements_to_plot):
            raise ValueError(f"Um ou mais IDs em ELEMENT_ID estão fora do intervalo válido [0, {max_id}].")
        
        # Define o título e a legenda para os elementos específicos
        title = f"Visualização do(s) Elemento(s) Específico(s): {elements_to_plot}"
        legend_patches = []
        if any(k in boundary_elements for k in elements_to_plot):
            legend_patches.append(Patch(facecolor='cyan', label='Elemento de Fronteira'))
        if any(k in internal_elements for k in elements_to_plot):
            legend_patches.append(Patch(facecolor='lightgreen', label='Elemento Interno'))
    else:
        # Comportamento original: usa plot_type para selecionar os elementos
        elements_to_plot, title, legend_patches = select_elements_to_plot(
            plot_type, boundary_elements, internal_elements)

    # 4. Procede com a plotagem se houver elementos a serem exibidos
    if not elements_to_plot:
        if ELEMENT_ID:
            print(f"IDs de elemento {ELEMENT_ID} não encontrados.")
        else:
            print(f"Nenhum elemento do tipo '{plot_type}' encontrado para plotar.")
        return

    print(f"\nPlotando {len(elements_to_plot)} elemento(s) (layout: '{layout}')...")

    # 5. Preparação de dados e configurações gerais para a plotagem
    plot_config = {
        'all_mesh_edges': extract_all_mesh_edges(dg.mesh.EToV),
        'domain_limits': (
            (dg.mesh.vx.min() - 0.1, dg.mesh.vx.max() + 0.1),
            (dg.mesh.vy.min() - 0.1, dg.mesh.vy.max() + 0.1),
            (dg.mesh.vz.min() - 0.1, dg.mesh.vz.max() + 0.1)
        ),
        'boundary_elements_set': boundary_elements,
        'boundary_node_ids': boundary_nodes,
        'layout': layout,
        'should_plot_nodes': not (layout == 'combined' and plot_type == 'all' and ELEMENT_ID is None),
        'alpha': 0.5 if layout == 'combined' else 0.3
    }
    
    # 6. Configuração da figura
    fig, axes = setup_figure(layout, len(elements_to_plot), title)
    
    # 7. Orquestração da plotagem
    if layout == 'separated':
        for i, k in enumerate(elements_to_plot):
            draw_single_element(axes[i], dg, k, plot_config)
    else: # 'combined'
        ax = axes[0]
        plot_wireframe(ax, dg.mesh.vx, dg.mesh.vy, dg.mesh.vz, plot_config['all_mesh_edges'])
        for k in elements_to_plot:
            # Corrigido para passar o argumento 'dg' que faltava na chamada anterior
            draw_single_element(ax, dg, k, plot_config)

    # 8. Finalização do plot
    if layout == 'separated':
        for j in range(len(elements_to_plot), len(axes)):
            axes[j].axis('off')
            
    if legend_patches:
        fig.legend(handles=legend_patches, loc='upper right', fontsize=12)
    plt.tight_layout(rect=[0, 0, 0.9, 0.95])


def plot_interface_elements(dg, group_id, partner_group_id, layout='combined'):
    """
    Plota os elementos de um grupo que participam de uma interface,
    destacando os nós de DG que estão exatamente na face da interface.
    A identificação dos elementos e nós é feita via mapas pré-calculados
    no objeto 'dg' (solver), e os nomes dos grupos são obtidos do objeto de malha.

    Parâmetros:
        dg: Objeto da classe Maxwell3D.
        group_id (int): O ID do grupo cujos elementos de interface serão plotados.
        partner_group_id (int): O ID do grupo vizinho na interface (para contexto no título).
        layout (str): 'separated' para um subplot por elemento ou 'combined'
                      para todos os elementos no mesmo subplot.
    """
    # 1. Obter nomes descritivos dos grupos a partir do objeto de malha
    group_name = dg.mesh.group_names.get(group_id, f"ID {group_id}")
    partner_name = dg.mesh.group_names.get(partner_group_id, f"ID {partner_group_id}")

    # 2. Obter dados da interface diretamente do objeto solver 'dg'
    try:
        elements_to_plot = getattr(dg, f'elements_G{group_id}')
        # Nós deste grupo que estão na interface
        nodes_to_highlight = set(getattr(dg, f'vmapIM_G{group_id}'))
    except AttributeError:
        print(f"\nErro: Mapas de interface para o grupo {group_id} ('{group_name}') não foram encontrados no objeto 'dg'.")
        print("-> Certifique-se de que 'dg.buildGroupInterfaceMaps()' foi executado corretamente.")
        print(f"-> Verifique também se uma interface entre os grupos {group_id} e {partner_group_id} foi descoberta em 'dg.interface_pairs'.")
        return

    if len(elements_to_plot) == 0:
        print(f"\nNenhum elemento de interface encontrado para o grupo '{group_name}'.")
        return

    print(f"\nPlotando {len(elements_to_plot)} elementos do Grupo '{group_name}' na interface com o Grupo '{partner_name}'...")

    # 3. Configurar a figura e os eixos com o título descritivo
    title = f"Elementos do Grupo '{group_name}' na Interface"
    fig, axes = setup_figure(layout, len(elements_to_plot), title)
    
    # 4. Preparar dados de fundo (wireframe da malha inteira)
    plot_config = {
        'all_mesh_edges': extract_all_mesh_edges(dg.mesh.EToV),
        'domain_limits': (
            (dg.mesh.vx.min() - 0.1, dg.mesh.vx.max() + 0.1),
            (dg.mesh.vy.min() - 0.1, dg.mesh.vy.max() + 0.1),
            (dg.mesh.vz.min() - 0.1, dg.mesh.vz.max() + 0.1)
        ),
        'layout': layout,
        'nodes_to_highlight': nodes_to_highlight, 
        'alpha': 0.5 if layout == 'combined' else 0.3
    }

    # 5. Orquestrar a plotagem dos elementos
    for i, k in enumerate(elements_to_plot):
        ax = axes[i] if layout == 'separated' else axes[0]
        
        vert_indices = dg.mesh.EToV[k]
        element_vertices = np.array([dg.mesh.vx[vert_indices], dg.mesh.vy[vert_indices], dg.mesh.vz[vert_indices]]).T
        node_coords_k = np.array([dg.x[:, k], dg.y[:, k], dg.z[:, k]]).T
        node_ids_k = dg.node_ids[:, k]

        plot_wireframe(ax, dg.mesh.vx, dg.mesh.vy, dg.mesh.vz, plot_config['all_mesh_edges'])
        plot_element_faces(ax, element_vertices, face_color='gold', alpha=plot_config['alpha'])
        plot_dg_nodes(ax, node_coords_k, node_ids_k, highlight_nodes=plot_config['nodes_to_highlight'])
        plot_element_label(ax, element_vertices, k)
        setup_plot_axes(ax, plot_config['domain_limits'])
        if layout == 'separated':
            ax.set_title(f"Elemento {k}", fontsize=10)

    # 6. Limpar e finalizar o plot
    if layout == 'separated':
        for j in range(len(elements_to_plot), len(axes)):
            axes[j].axis('off')

    # Adiciona uma legenda mais descritiva
    legend_patch = Patch(facecolor='magenta', edgecolor='black', label=f"Nós do Grupo '{group_name}' na Interface")
    fig.legend(handles=[legend_patch], loc='upper right')
    plt.tight_layout(rect=[0, 0, 0.9, 0.95])


def plot_interface_nodes(dg, title="Visualização dos Nós da Interface"):
    """
    Gera um gráfico 3D dos nós que compõem as interfaces físicas.
    Os pares de IDs de nós (Minus, Plus) são listados em uma caixa de
    legenda externa para evitar sobreposição.

    Parâmetros
    ----------
    dg : Maxwell3D
        A instância da classe de discretização espacial que contém os nós e mapas.
    title : str, opcional
        O título do gráfico.
    """
    fig = plt.figure(figsize=(16, 10))
    ax = fig.add_subplot(111, projection='3d')

    # 1. Plota todos os nós da malha como um fundo cinza para dar contexto
    ax.scatter(dg.x, dg.y, dg.z, c='gray', s=2, alpha=0.1, label='Todos os Nós da Malha')

    # 2. Verifica se os pares de interface foram calculados
    if not hasattr(dg, 'interface_pairs') or not dg.interface_pairs:
        print("Nenhuma interface física foi descoberta. Nenhum mapa de grupo foi construído.")
        ax.set_title("Nenhuma interface encontrada", fontsize=16)
        plt.tight_layout()
        return
    
    # --- Pré-calcular o total de nós ---
    total_nodes_plotted = 0
    for id_A, id_B in dg.interface_pairs:
        try:
            nodes_M = getattr(dg, f'vmapIM_G{id_A}')
            total_nodes_plotted += len(nodes_M)
        except AttributeError:
            continue # Pula se o mapa não for encontrado para este par

    # 3. Itera sobre cada par de interface para plotar os nós e coletar os rótulos
    colors = plt.cm.get_cmap('viridis', len(dg.interface_pairs))
    
    # Monta o título da legenda com a contagem total
    legend_title = f"Mapeamento de {total_nodes_plotted} Nós (Minus, Plus)"
    legend_text_parts = [legend_title, "------------------------------------"]

    print("\nPlotando nós de interface e preparando a legenda...")
    for i, (id_A, id_B) in enumerate(dg.interface_pairs):
        try:
            nodes_M = getattr(dg, f'vmapIM_G{id_A}')
            nodes_P = getattr(dg, f'vmapIP_G{id_A}')
        except AttributeError:
            print(f"⚠️ Aviso: Não foi possível encontrar mapas de interface para o par ({id_A}, {id_B}). Pulando.")
            continue

        if len(nodes_M) == 0:
            continue

        # Adiciona o cabeçalho da interface atual à lista de textos da legenda
        group_A_name = dg.mesh.group_names.get(id_A, f'Grupo {id_A}')
        group_B_name = dg.mesh.group_names.get(id_B, f'Grupo {id_B}')
        legend_text_parts.append(f"\nInterface '{group_A_name}' <> '{group_B_name}':")

        # Coleta as coordenadas e plota os pontos
        x_coords = dg.x.ravel('F')[nodes_M]
        y_coords = dg.y.ravel('F')[nodes_M]
        z_coords = dg.z.ravel('F')[nodes_M]
        color = colors(i / len(dg.interface_pairs) if len(dg.interface_pairs) > 1 else 0.5)
        ax.scatter(x_coords, y_coords, z_coords, c=[color], s=40, depthshade=True, edgecolor='k', linewidth=0.5)

        # Formata os pares de IDs para a legenda externa
        current_pairs = [f"({m}, {p})" for m, p in zip(nodes_M, nodes_P)]
        
        # Organiza os pares em linhas para melhor legibilidade
        pairs_per_line = 4
        for chunk_idx in range(0, len(current_pairs), pairs_per_line):
            chunk = current_pairs[chunk_idx:chunk_idx + pairs_per_line]
            line = ", ".join(chunk)
            legend_text_parts.append("  " + line) # Indentação para clareza

        print(f"🔎 Interface entre Grupo {id_A} e {id_B}: {len(nodes_M)} nós plotados.")

    # 4. Configurações finais do gráfico
    ax.set_title(title, fontsize=16)
    ax.set_xlabel("Eixo X")
    ax.set_ylabel("Eixo Y")
    ax.set_zlabel("Eixo Z")
    ax.view_init(elev=30, azim=-60)

    # 5. Adiciona a caixa de texto com todos os mapeamentos de ID
    final_legend_text = "\n".join(legend_text_parts)
    fig.text(0.72, 0.9, final_legend_text,
             transform=fig.transFigure,
             fontsize=9,
             family='monospace',
             va='top',
             ha='left',
             bbox=dict(boxstyle='round,pad=0.5', fc='aliceblue', ec='black', alpha=0.95))

    # Ajusta o layout para criar espaço para a caixa de legenda
    plt.tight_layout(rect=[0, 0, 0.7, 1])

