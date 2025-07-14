import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
from matplotlib.collections import LineCollection

from .mesh2d import Mesh2D
from .dg2d import Maxwell2D

# NOVA CLASSE FILHA: Responsável pela visualização
class Mesh2DVisualizer(Mesh2D):
    def __init__(self, vx, vy, EToV, physical_groups=None, boundary_label="PEC"):
        # Chama o construtor da classe pai (Mesh2D) para inicializar os dados
        super().__init__(vx, vy, EToV, physical_groups=physical_groups, boundary_label=boundary_label)

        # Atributos de estilo para customização dos plots
        self.figure_size = (8, 8)
        self.title_font_size = 10
        self.element_label_font_size = 10
        self.node_label_font_size = 7
        self.node_marker_color = 'red'
        self.element_label_color = 'darkslateblue'

    def vertices_and_elements(self, show_vertices=True, title='Malha Triangular'):
        """ Plota a malha triangular com numeração dos elementos e vértices. """
        # Esta classe herda Mesh2D
        triangulation = self.getTriangulation()
        
        fig, ax = plt.subplots(figsize=self.figure_size)
        ax.set_aspect('equal')
        ax.triplot(triangulation, 'k-', lw=1, label='Elementos')

        if show_vertices:
            ax.plot(self.vx, self.vy, 'o', markersize=4, label='Vértices')

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
        ax.legend(loc='upper right', bbox_to_anchor=(1, 1.1))
        plt.grid(True, linestyle='--', alpha=0.5)


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
        ax.legend(handles=legend_handles, labels=legend_labels, loc='upper right', bbox_to_anchor=(1, 1.1))
        plt.grid(True, linestyle='--', alpha=0.5)


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
    
    
    def _plot_collocation_data(self, ax, dg, map_type, shift_distance, marker_label):
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
            ax.plot(marker_x_coords, marker_y_coords, 'o', markersize=5, color=self.node_marker_color, alpha=0.8, linestyle='none', label=marker_label)


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
        
        ax.triplot(self.getTriangulation(), 'k-', lw=0.5, alpha=0.3, label='_nolegend_')

        dynamic_title = title.replace('(vmap)', f'(vmap{map_type})')
        marker_label = f'Nós de Colocação ({map_type})'

        # Chama os métodos auxiliares para realizar a plotagem
        self._plot_element_indices(ax)
        self._plot_collocation_data(ax, dg, map_type, shift_distance, marker_label)

        # Finaliza a configuração do gráfico
        ax.set_title(dynamic_title, fontsize=self.title_font_size, fontweight='bold')
        ax.set_xlabel('x (m)')
        ax.set_ylabel('y (m)')
        ax.legend(loc='upper right', bbox_to_anchor=(1, 1.1))
        plt.grid(True, linestyle='--', alpha=0.3)

