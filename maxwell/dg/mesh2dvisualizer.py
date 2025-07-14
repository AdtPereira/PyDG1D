import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
from matplotlib.collections import LineCollection

from .mesh2d import Mesh2D


# NOVA CLASSE FILHA: Responsável pela visualização
class Mesh2DVisualizer(Mesh2D):
    def __init__(self, vx, vy, EToV, boundary_label="PEC"):
        # Chama o construtor da classe pai (Mesh2D) para inicializar os dados
        super().__init__(vx, vy, EToV, boundary_label)


    def plot_vertices_and_elements(self, show_vertices=True, title='Malha Triangular'):
        """ Plota a malha triangular com numeração dos elementos e vértices. """
        # Esta classe herda Mesh2D
        triangulation = self.getTriangulation()
        
        fig, ax = plt.subplots()
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

        ax.set_title(title, fontsize=10)
        ax.set_xlabel('x (m)')
        ax.set_ylabel('y (m)')
        ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
        plt.grid(True, linestyle='--', alpha=0.5)


    def plot_physical_group_map(self, physical_group_data: dict, title='Mapa de Grupos Físicos'):
        """
        Plota a malha, tratando e colorindo os elementos de cada grupo físico
        de acordo com sua dimensão (1D para linhas, 2D para superfícies).
        """
        fig, ax = plt.subplots()
        ax.set_aspect('equal')
        
        colors = plt.cm.get_cmap('viridis', len(physical_group_data))

        # Listas para criar a legenda manualmente no final
        legend_handles = []
        legend_labels = []

        for i, (group_name, data) in enumerate(physical_group_data.items()):
            group_dim = data['dim']
            vx_group = data['vx']
            vy_group = data['vy']
            EToV = data['EToV']

            # Converte a conectividade para índices locais (base 0)
            unique_nodes, inverse_indices = np.unique(EToV.flatten(), return_inverse=True)
            connectivity_local = inverse_indices.reshape(EToV.shape)

            # Adiciona o nome do grupo e a cor à legenda
            legend_labels.append(group_name)
            
            # --- Lógica de Plotagem baseada na Dimensão ---
            if group_dim == 2:
                # Para elementos 2D, usa Triangulation e tripcolor
                triangulation = mtri.Triangulation(vx_group, vy_group, connectivity_local)
                ax.tripcolor(triangulation, facecolors=np.repeat(i, len(connectivity_local)),
                            cmap='viridis', vmin=0, vmax=len(physical_group_data)-1,
                            alpha=0.6)
                ax.triplot(triangulation, 'k-', lw=0.5)
                legend_handles.append(plt.Rectangle((0,0),1,1, color=colors.colors[i], alpha=0.6))

            elif group_dim == 1:
                # Para elementos 1D, usa LineCollection
                # Cria uma lista de segmentos de linha: [ ((x1, y1), (x2, y2)), ... ]
                segments = [((vx_group[p1], vy_group[p1]), (vx_group[p2], vy_group[p2])) for p1, p2 in connectivity_local]
                lines = LineCollection(segments, colors=colors.colors[i], linewidths=2.5, label=group_name)
                ax.add_collection(lines)
                legend_handles.append(lines)

        ax.set_title(title, fontsize=10, fontweight='bold')
        ax.set_xlabel('x (m)')
        ax.set_ylabel('y (m)')
        ax.legend(handles=legend_handles, labels=legend_labels, loc='center left', bbox_to_anchor=(1, 0.5))
        plt.grid(True, linestyle='--', alpha=0.5)