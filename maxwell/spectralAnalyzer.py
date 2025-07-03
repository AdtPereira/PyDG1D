import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm

class ResonantCavity3D:
    def __init__(self, problem):
        self.DIM = problem['DIM']
        self.m = problem['m']
        self.n = problem['n']
        self.L = problem['L']
        self.t0 = problem['t0']
        
        self.t_final = 10
        self.kx = self.m * np.pi / self.L
        self.ky = self.n * np.pi / self.L
        self.w = np.sqrt(self.kx**2 + self.ky**2)

    def Ez_field(self, x, y, z, t):
        return np.sin(self.kx * x) * np.sin(self.ky * y) * np.cos(self.w * t)


    def Hx_field(self, x, y, z, t):
        return - (self.ky / self.w) * np.sin(self.kx * x) * np.cos(self.ky * y) * np.sin(self.w * t)


    def Hy_field(self, x, y, z, t):
        return (self.kx / self.w) * np.cos(self.kx * x) * np.sin(self.ky * y) * np.sin(self.w * t)


    def dyEz_field(self, x, y, z, t):
        return (self.ky) * np.sin(self.kx * x) * np.cos(self.ky * y) * np.cos(self.w * t)


    def dxEz_field(self, x, y, z, t):
        return (self.kx) * np.cos(self.kx * x) * np.sin(self.ky * y) * np.cos(self.w * t)


    # def plot_dg_solution(self, mesh, x_i, y_i, z_i, ez_i, title, slice_tol=0.05):
    #     """
    #     Plota a solução 3D interpolada e um corte 2D no plano z=0.
    #     """
    #     fig = plt.figure(figsize=(20, 9))
    #     fig.suptitle(title, fontsize=18)

    #     # --- Subplot 1: Visualização 3D completa ---
    #     ax1 = fig.add_subplot(1, 2, 1, projection='3d')
    #     scatter = ax1.scatter(x_i, y_i, z_i, c=ez_i, cmap='viridis', s=5, alpha=0.6)
    #     cbar = fig.colorbar(scatter, ax=ax1, shrink=0.7, aspect=15, label="Amplitude do Campo (Ez)")
        
    #     ax1.set_title("Visualização 3D Completa")
    #     ax1.set_xlabel('X')
    #     ax1.set_ylabel('Y')
    #     ax1.set_zlabel('Z')
        
    #     # --- Subplot 2: Corte no plano z=0 ---
    #     ax2 = fig.add_subplot(1, 2, 2)
        
    #     # Filtra os pontos que estão próximos ao plano z=0
    #     mask = np.abs(z_i) < slice_tol
    #     x_slice, y_slice, ez_slice = x_i[mask], y_i[mask], ez_i[mask]
        
    #     # Verifica se há pontos suficientes para plotar
    #     if len(x_slice) > 3:
    #         # tricontourf é excelente para dados não estruturados
    #         contour = ax2.tricontourf(x_slice, y_slice, ez_slice, levels=20, cmap='viridis')
    #         fig.colorbar(contour, ax=ax2, label="Amplitude do Campo (Ez)")
    #         ax2.set_title(f"Corte da Solução no Plano z ≈ 0 (±{slice_tol})")
    #     else:
    #         ax2.text(0.5, 0.5, 'Nenhum dado encontrado para o corte.', ha='center', va='center')
    #         ax2.set_title("Corte da Solução no Plano z ≈ 0")
            
    #     ax2.set_xlabel('X')
    #     ax2.set_ylabel('Y')
    #     ax2.set_aspect('equal', adjustable='box')
    #     ax2.grid(True, linestyle='--', alpha=0.5)

    #     plt.tight_layout(rect=[0, 0, 1, 0.96])


    # def plot_analytical_solution(self, mesh, x_interp, y_interp, z_interp, t, title):
        """
        Plota a solução analítica em um grid de pontos fornecido, usando o mesmo
        estilo da plot_interpolated_solution para permitir uma comparação direta.

        Parâmetros
        ----------
        mesh : Mesh3D
            O objeto de malha, usado para desenhar as arestas dos elementos.
        x_interp, y_interp, z_interp : ndarray
            Coordenadas dos pontos do grid de interpolação.
        t : float
            O instante de tempo no qual a solução analítica será calculada.
        title : str
            Título do gráfico.
        """
        fig = plt.figure(figsize=(12, 9))
        ax = fig.add_subplot(111, projection='3d')

        # 1. Calcular o campo analítico nos pontos do grid de interpolação
        ez_analytical = self.Ez_field(x_interp, y_interp, z_interp, t)

        # 2. Plotar a solução analítica como uma nuvem de pontos colorida
        scatter = ax.scatter(x_interp, y_interp, z_interp, c=ez_analytical, cmap='viridis', s=5, alpha=0.7)
        
        # 3. Adicionar uma barra de cores
        cbar = fig.colorbar(scatter, ax=ax, shrink=0.7, aspect=10)
        cbar.set_label("Amplitude Analítica do Campo (Ez)")

        # 4. Sobrepor as arestas da malha para dar contexto (código reutilizado)
        for k in range(mesh.number_of_elements()):
            vert_indices = mesh.EToV[k, :]
            verts = np.array([mesh.vx[vert_indices], mesh.vy[vert_indices], mesh.vz[vert_indices]]).T
            
            # Arestas de um tetraedro
            edges = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]
            for i, j in edges:
                ax.plot([verts[i, 0], verts[j, 0]], 
                        [verts[i, 1], verts[j, 1]], 
                        [verts[i, 2], verts[j, 2]], 'k-', lw=0.5, alpha=0.3)

        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.set_zlabel('Z')
        ax.set_title(title, fontsize=16)
        ax.view_init(elev=25, azim=-135) # Ajusta o ângulo da câmera


    def plot_field(self, x_i, y_i, z_i, field_data, title, slice_tol=0.05, cbar_label="Amplitude do Campo (Ez)"):
        """
        Plota um campo escalar 3D e um corte 2D no plano z=0.

        Esta função é genérica e pode ser usada tanto para a solução numérica
        interpolada quanto para a solução analítica calculada nos mesmos pontos.

        Parâmetros
        ----------
        x_i, y_i, z_i : ndarray
            Coordenadas dos pontos do grid.
        field_data : ndarray
            Valores do campo escalar a serem plotados (ex: Ez interpolado ou analítico).
        title : str
            Título geral para a figura.
        slice_tol : float, opcional
            Tolerância para selecionar pontos para o corte no plano z=0.
        cbar_label : str, opcional
            Rótulo para a barra de cores.
        """
        fig = plt.figure(figsize=(20, 9))
        fig.suptitle(title, fontsize=18)

        # --- Subplot 1: Visualização 3D completa ---
        ax1 = fig.add_subplot(1, 2, 1, projection='3d')
        # Usa 'field_data' para a cor dos pontos
        scatter = ax1.scatter(x_i, y_i, z_i, c=field_data, cmap='viridis', s=5, alpha=0.6)
        cbar = fig.colorbar(scatter, ax=ax1, shrink=0.7, aspect=15, label=cbar_label)

        ax1.set_title("Visualização 3D Completa")
        ax1.set_xlabel('X')
        ax1.set_ylabel('Y')
        ax1.set_zlabel('Z')
        ax1.view_init(elev=25, azim=-135)

        # --- Subplot 2: Corte no plano z=0 ---
        ax2 = fig.add_subplot(1, 2, 2)

        # Filtra os pontos que estão próximos ao plano z=0
        mask = np.abs(z_i) < slice_tol
        x_slice, y_slice, field_slice = x_i[mask], y_i[mask], field_data[mask]

        # Verifica se há pontos suficientes para plotar
        if len(x_slice) > 3:
            # tricontourf é excelente para dados não estruturados
            contour = ax2.tricontourf(x_slice, y_slice, field_slice, levels=20, cmap='viridis')
            fig.colorbar(contour, ax=ax2, label=cbar_label)
            ax2.set_title(f"Corte da Solução no Plano z ≈ 0 (±{slice_tol})")
        else:
            ax2.text(0.5, 0.5, 'Nenhum dado encontrado para o corte.', ha='center', va='center')
            ax2.set_title("Corte da Solução no Plano z ≈ 0")

        ax2.set_xlabel('X')
        ax2.set_ylabel('Y')
        ax2.set_aspect('equal', adjustable='box')
        ax2.grid(True, linestyle='--', alpha=0.5)

        plt.tight_layout(rect=[0, 0, 1, 0.96])