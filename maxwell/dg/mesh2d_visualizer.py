import numpy as np
import os
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
from matplotlib.collections import LineCollection

from .mesh2d import Mesh2D
from .dg2d import Maxwell2D

class Mesh2DVisualizer(Mesh2D):
    def __init__(self, vx, vy, EToV, physical_groups=None, boundary_label="PEC", filePath=None):
        # Chama o construtor da classe pai (Mesh2D) para inicializar os dados
        super().__init__(vx, vy, EToV, physical_groups=physical_groups, boundary_label=boundary_label)

        # Atributos de estilo para customização dos plots
        self.filePath = filePath                    # Caminho para salvar o gráfico, se necessário
        self.bbox_anchor = (1.0, 1.15)              # Posição da legenda fora do gráfico
        self.figure_size = (8, 8)                   # Tamanho da figura do gráfico
        self.title_font_size = 10                   # Tamanho da fonte do título do gráfico
        self.element_label_font_size = 10           # Tamanho da fonte dos rótulos dos elementos
        self.node_label_font_size = 7               # Tamanho da fonte dos rótulos dos nós
        self.marker_size = 3                        # Tamanho dos marcadores dos nós
        self.node_marker_color = 'red'              # Cor dos marcadores dos nós
        self.element_label_color = 'darkslateblue'  # Cor dos rótulos dos elementos


    def vertices_and_elements(self, show_vertices=True, title='Malha Triangular'):
        """ Plota a malha triangular com numeração dos elementos e vértices. """
        # Esta classe herda Mesh2D
        triangulation = self.getTriangulation()
        
        fig, ax = plt.subplots(figsize=self.figure_size)
        ax.set_aspect('equal')
        ax.triplot(triangulation, 'k-', lw=1, label='Elementos')

        if show_vertices:
            ax.plot(self.vx, self.vy, 'o', markersize=self.marker_size, label='Vértices')

            # Itera sobre o número total de vértices
            for i in range(self.number_of_vertices()):
                ax.text(self.vx[i], self.vy[i], f' {i}',# O texto a ser exibido (o índice do vértice)
                        fontsize=9,                     # Tamanho da fonte para não poluir o gráfico
                        color='blue',                   # Cor do texto para diferenciar dos vértices
                        fontweight='bold',              # Deixa o texto mais destacado
                        verticalalignment='bottom',     # Alinha o texto acima do ponto
                        horizontalalignment='left')     # Alinha o texto à direita do ponto
                
            # Itera sobre cada elemento
            for k in range(self.number_of_elements()):
                ax.text(np.mean(self.vx[self.EToV[k]]),
                        np.mean(self.vy[self.EToV[k]]),
                        str(k),
                        fontsize=9,
                        color='red',
                        fontweight='bold',
                        horizontalalignment='center',
                        verticalalignment='center')

        ax.set_title(title, fontsize=self.title_font_size)
        ax.set_xlabel('x (m)')
        ax.set_ylabel('y (m)')
        ax.legend(loc='upper right', bbox_to_anchor=self.bbox_anchor)
        plt.grid(True, linestyle='--', alpha=0.5)

        full_path = self.filePath + "_vertices.svg"
        plt.savefig(full_path, dpi=300, bbox_inches='tight', pad_inches=0.1)
        print(f"\n📄 Graphics saved at {full_path}")


    def gmsh_physical_group_map(self, physical_group_data: dict, title='Mapa de Grupos Físicos'):
        """
        Plota a malha, tratando e colorindo os elementos de cada grupo físico
        de acordo com sua dimensão (1D para linhas, 2D para superfícies).
        """
        fig, ax = plt.subplots(figsize=self.figure_size)
        ax.set_aspect('equal')
        
        colors = plt.cm.get_cmap('viridis', len(physical_group_data))

        # Listas para criar a legenda manualmente no final
        legend_handles = []
        legend_labels = []

        for i, (TagGroup, DataGroup) in enumerate(physical_group_data.items()):
            dim = DataGroup['dim']
            vx = DataGroup['vx']
            vy = DataGroup['vy']
            EToV = DataGroup['EToV']
            name = DataGroup['name']

            # Converte a conectividade para índices locais (base 0)
            unique_nodes, inverse_indices = np.unique(EToV.flatten(), return_inverse=True)
            connectivity_local = inverse_indices.reshape(EToV.shape)

            # Adiciona o nome do grupo e a cor à legenda
            legend_labels.append(f"{name} ({TagGroup})")

            # --- Lógica de Plotagem baseada na Dimensão ---
            if dim == 2:
                # Para elementos 2D, usa Triangulation e tripcolor
                triangulation = mtri.Triangulation(vx, vy, connectivity_local)
                ax.tripcolor(triangulation, facecolors=np.repeat(i, len(connectivity_local)),
                            cmap='viridis', vmin=0, vmax=len(physical_group_data)-1,
                            alpha=0.6)
                ax.triplot(triangulation, 'k-', lw=0.5)
                legend_handles.append(plt.Rectangle((0,0),1,1, color=colors.colors[i], alpha=0.6))

            elif dim == 1:
                # Para elementos 1D, usa LineCollection
                # Cria uma lista de segmentos de linha: [ ((x1, y1), (x2, y2)), ... ]
                segments = [((vx[p1], vy[p1]), (vx[p2], vy[p2])) for p1, p2 in connectivity_local]
                lines = LineCollection(segments, colors=colors.colors[i], linewidths=2.5, label=TagGroup)
                ax.add_collection(lines)
                legend_handles.append(lines)

        ax.set_title(title, fontsize=self.title_font_size)
        ax.set_xlabel('x (m)')
        ax.set_ylabel('y (m)')
        ax.legend(handles=legend_handles, labels=legend_labels, loc='upper right', bbox_to_anchor=self.bbox_anchor)
        plt.grid(True, linestyle='--', alpha=0.5)

        full_path = self.filePath + "_gmsh_groups.svg"
        plt.savefig(full_path, dpi=300, bbox_inches='tight', pad_inches=0.1)
        print(f"\n📄 Graphics saved at {full_path}")


    def _plot_element_indices(self, ax):
        """Método auxiliar para plotar os índices 'k' no centro de cada elemento."""
        centroids_x = np.array([np.mean(self.vx[self.EToV[k]]) for k in range(self.number_of_elements())])
        centroids_y = np.array([np.mean(self.vy[self.EToV[k]]) for k in range(self.number_of_elements())])

        for k in range(self.number_of_elements()):
            ax.text(centroids_x[k], centroids_y[k], str(k),
                    fontsize=self.element_label_font_size,
                    color=self.element_label_color,
                    fontweight='bold', ha='center', va='center', alpha=0.6)
        return centroids_x, centroids_y
    
    
    def _plot_interface_dofs(self, ax, dg, map_type, shift_distance, marker_label):
        """Método auxiliar para plotar os nós e os índices dos pontos de colocação."""
        vmap_to_plot = dg.vmapM if map_type == 'M' else dg.vmapP
        vmap_reshaped = vmap_to_plot.reshape((dg.n_fp, dg.n_faces, self.number_of_elements()), order='F')
        
        x_coords_flat = dg.x.ravel('F')
        y_coords_flat = dg.y.ravel('F')
        EToE, _ = self.connectivityMatrices()
        
        centroids_x = np.array([np.mean(self.vx[self.EToV[k]]) for k in range(self.number_of_elements())])
        centroids_y = np.array([np.mean(self.vy[self.EToV[k]]) for k in range(self.number_of_elements())])

        marker_x_coords, marker_y_coords = [], []

        for k in range(self.number_of_elements()):
            for f in range(dg.n_faces):
                owner_element_idx = k if map_type == 'M' else EToE[k, f]
                centroid_x, centroid_y = centroids_x[owner_element_idx], centroids_y[owner_element_idx]

                for p in range(dg.n_fp):
                    global_idx = vmap_reshaped[p, f, k]
                    px, py = x_coords_flat[global_idx], y_coords_flat[global_idx]
                    marker_x_coords.append(px)
                    marker_y_coords.append(py)
                    
                    vec_x, vec_y = centroid_x - px, centroid_y - py
                    norm = np.sqrt(vec_x**2 + vec_y**2) + 1e-9
                    text_x, text_y = px + shift_distance * (vec_x / norm), py + shift_distance * (vec_y / norm)
                    
                    ax.text(text_x, text_y, str(global_idx), fontsize=self.node_label_font_size, color='black', fontweight='bold', ha='center', va='center')
        
        if marker_x_coords:
            ax.plot(marker_x_coords, marker_y_coords, 'o', markersize=self.marker_size, color=self.node_marker_color, alpha=0.8, linestyle='none', label=marker_label)


    def _plot_internal_dofs(self, ax, dg, shift_distance):
        """Método auxiliar para plotar os nós e os índices dos pontos de colocação INTERNOS."""
        x_coords_flat = dg.x.ravel('F')
        y_coords_flat = dg.y.ravel('F')
        
        # Identifica os índices dos nós que não estão nas faces (vmapM)
        face_indices = set(dg.vmapM.flatten())
        all_indices = set(range(len(x_coords_flat)))
        internal_indices = sorted(list(all_indices - face_indices))
        
        if not internal_indices:
            return # Não há pontos internos para plotar

        # Calcula os centroides para o deslocamento do texto
        centroids_x = np.array([np.mean(self.vx[self.EToV[k]]) for k in range(self.number_of_elements())])
        centroids_y = np.array([np.mean(self.vy[self.EToV[k]]) for k in range(self.number_of_elements())])
        
        marker_x_coords, marker_y_coords = [], []

        for global_idx in internal_indices:
            # Encontra a qual elemento 'k' o ponto pertence
            k = global_idx // dg.number_of_nodes_per_element()
            
            px, py = x_coords_flat[global_idx], y_coords_flat[global_idx]
            marker_x_coords.append(px)
            marker_y_coords.append(py)
            
            # Calcula o vetor do ponto ao centroide do elemento para deslocar o texto
            centroid_x, centroid_y = centroids_x[k], centroids_y[k]
            vec_x, vec_y = centroid_x - px, centroid_y - py
            norm = np.sqrt(vec_x**2 + vec_y**2) + 1e-9
            text_x, text_y = px + shift_distance * (vec_x / norm), py + shift_distance * (vec_y / norm)
            
            ax.text(text_x, text_y, str(global_idx), fontsize=self.node_label_font_size, color='darkgreen', ha='center', va='center')
        
        ax.plot(marker_x_coords, marker_y_coords, 's', markersize=self.marker_size, color='green', alpha=0.8, linestyle='none', label='Nós Internos')


    def _plot_physical_group_edges(self, ax):
        """
        Método auxiliar para sobrepor as arestas que pertencem a grupos físicos 1D.

        Args:
            ax (matplotlib.axes.Axes): O eixo onde a plotagem será realizada.
        """
        # Verifica se existem grupos físicos para evitar erros
        if not self.physical_groups:
            return
        
        colors = plt.cm.get_cmap('viridis', len(self.physical_groups))

        # Itera sobre todos os grupos físicos definidos
        for i, (key, DataGroup) in enumerate(self.physical_groups.items()):
            # Plota apenas os grupos que são 1D (linhas/arestas)
            if DataGroup.get('dim') == 1:
                vx_group = DataGroup['vx']
                vy_group = DataGroup['vy']
                EToV_group = DataGroup['EToV']
                name_group = f"{DataGroup['name']} (Tag {key})"

                # Lida com a conectividade local do grupo
                unique_nodes, inverse_indices = np.unique(EToV_group.flatten(), return_inverse=True)
                connectivity_local = inverse_indices.reshape(EToV_group.shape)

                # Cria uma lista de segmentos de linha a partir da conectividade
                # O formato é [((x1, y1), (x2, y2)), ...]
                segments = [((vx_group[p1], vy_group[p1]), (vx_group[p2], vy_group[p2])) for p1, p2 in connectivity_local]

                # Usa LineCollection para plotar as arestas de forma eficiente
                # com uma cor e espessura distintas para destaque
                lines = LineCollection(segments, colors=colors(i), linewidths=1,
                                       zorder=5, label=name_group)
                ax.add_collection(lines)


    def collocation_points(self, dg: Maxwell2D, map_type='M', shift_distance=0.08, title='Índices Globais dos Pontos de Colocação (vmap)'):
        """
        Orquestra a plotagem dos nós, índices de elementos e índices de pontos de colocação.

        Args:
            dg (Maxwell2D): O objeto do solver DG.
            map_type (str): O tipo de mapa a ser plotado, 'M' para vmapM ou 'P' para vmapP.
            shift_distance (float): Fator de deslocamento do texto para dentro do elemento.
            title (str): Título do gráfico.
        """
        if map_type not in ['M', 'P']:
            raise ValueError("O argumento map_type deve ser 'M' ou 'P'.")

        fig, ax = plt.subplots(figsize=self.figure_size)
        ax.set_aspect('equal')
        
        # Desenha a malha base
        ax.triplot(self.getTriangulation(), 'k-', lw=0.5, alpha=0.3, label='_nolegend_')

        # Sobrepõe as arestas dos grupos físicos para destaque
        self._plot_physical_group_edges(ax)

        # Plota os índices dos elementos
        dynamic_title = title.replace('(vmap)', f'(vmap{map_type})')
        marker_label = f'Nós de Colocação ({map_type})'

        # Chama os métodos auxiliares para realizar a plotagem
        self._plot_element_indices(ax)
        self._plot_interface_dofs(ax, dg, map_type, shift_distance, marker_label)
        self._plot_internal_dofs(ax, dg, shift_distance)

        # Finaliza a configuração do gráfico
        ax.set_title(dynamic_title, fontsize=self.title_font_size, fontweight='bold')
        ax.set_xlabel('x (m)')
        ax.set_ylabel('y (m)')
        ax.legend(loc='upper right', bbox_to_anchor=self.bbox_anchor)
        plt.grid(True, linestyle='--', alpha=0.3)

        full_path = self.filePath + "_collocation_points.svg"
        plt.savefig(full_path, dpi=300, bbox_inches='tight', pad_inches=0.1)
        print(f"\n📄 Graphics saved at {full_path}")

