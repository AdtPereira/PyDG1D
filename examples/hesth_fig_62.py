import numpy as np
import matplotlib.pyplot as plt
from matplotlib.tri import Triangulation
from scipy.special import eval_jacobi
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from mpl_toolkits.axes_grid1 import make_axes_locatable

# Define the orthonormal basis function from Eq. (6.6)
def psi_ij(r, s, i, j):
    a = 2 * (1 + r) / (1 - s + 1e-14) - 1
    b = s
    P_i = eval_jacobi(i, 0, 0, a)
    P_j = eval_jacobi(j, 2*i + 1, 0, b)
    return np.sqrt(2) * P_i * P_j * (1 - b)**i

# Define the number of points in each direction for plotting
Nplot = 40
r = np.linspace(-1, 1, Nplot)
s = np.linspace(-1, 1, Nplot)
R, S = np.meshgrid(r, s)

# Mask out points outside the reference triangle
mask = (R + S <= 0)
R = R[mask]
S = S[mask]

# Prepare triangulation
tri = Triangulation(R, S)

# Define maximum order
N_max = 2

# Define edges of the reference triangle
edges = np.array([
    [[-1, -1], [1, -1]],
    [[1, -1], [-1, 1]],
    [[-1, 1], [-1, -1]]
])

# Adiciona os pontos nodais e o índice da função base nas figuras
# Define os pontos nodais do triângulo de referência para visualização
def reference_triangle_nodes(N):
    L = []
    for n in range(N+1):
        for m in range(N+1 - n):
            r = -1 + 2 * m / N
            s = -1 + 2 * n / N
            L.append((r, s))
    return np.array(L)

for N in range(N_max + 1):
    fig = plt.figure(figsize=(3 * (N + 1), 6))
    fig.suptitle(f'Ordem N = {N}', fontsize=16)

    nodes = reference_triangle_nodes(N) if N > 0 else np.array([[ -1, -1 ]])

    for idx, i in enumerate(range(N + 1)):
        j = N - i
        Z = psi_ij(R, S, i, j)

        # Heatmap
        ax1 = fig.add_subplot(2, N + 1, idx + 1)
        tpc = ax1.tripcolor(R, S, Z, shading='gouraud', cmap='coolwarm')
        ax1.plot(nodes[:, 0], nodes[:, 1], 'ko', markersize=3)
        ax1.set_aspect('equal')
        ax1.set_xticks([])
        ax1.set_yticks([])
        ax1.set_title(f'Heatmap $(i,j)=({i},{j})$', fontsize=10)
        divider = make_axes_locatable(ax1)
        cax = divider.append_axes("right", size="5%", pad=0.05)
        fig.colorbar(tpc, cax=cax)

        # Surface
        ax2 = fig.add_subplot(2, N + 1, N + 1 + idx + 1, projection='3d')
        ax2.plot_trisurf(R, S, Z, triangles=tri.triangles, cmap='coolwarm', edgecolor='none')

        # Arestas do triângulo de referência
        for edge in edges:
            x_edge, y_edge = edge[:, 0], edge[:, 1]
            z_edge = np.zeros_like(x_edge)
            ax2.plot(x_edge, y_edge, z_edge, color='k', linewidth=1.5)

        # Pontos nodais
        ax2.scatter(nodes[:, 0], nodes[:, 1], np.zeros_like(nodes[:, 0]), color='k', s=10)

        # Índice da base
        ax2.text2D(0.05, 0.9, f"$\\psi_{{({i},{j})}}$", transform=ax2.transAxes, fontsize=10)

        ax2.set_xticks([])
        ax2.set_yticks([])
        ax2.set_zticks([])
        ax2.view_init(elev=30, azim=45)
        ax2.set_title(f'Superfície $(i,j)=({i},{j})$', fontsize=10)

    plt.tight_layout()
    plt.show()



