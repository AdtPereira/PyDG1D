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
            ax.plot(x, y, z, 'k--', linewidth=0.8, zorder=0)

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