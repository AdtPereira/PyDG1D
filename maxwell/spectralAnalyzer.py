import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation


class ResonantCavity3D:
    def __init__(self, problem):
        self.DIM = problem['DIM']
        self.m = problem['m']
        self.n = problem['n']
        self.Lx = problem['domain']['Lx']
        self.Ly = problem['domain']['Ly']
        self.Lz = problem['domain']['Lz']
        self.t0 = problem['dg']['t0']        
        self.t_final = problem['dg']['t_final']
        
        self.kx = self.m * np.pi / self.Lx
        self.ky = self.n * np.pi / self.Ly
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


    def plot_L2_error(self, error_data) -> None:
        """
        Plota a evolução do erro L2 ao longo do tempo para um ou mais campos.

        A função gera um gráfico com o tempo no eixo x e o erro L2 no eixo y,
        utilizando uma escala logarítmica para melhor visualização de pequenas
        variações de erro.

        Parâmetros
        ----------
        error_data : dict
            Dicionário contendo os dados do erro. Espera-se a seguinte estrutura:
            {
                'time': [t_0, t_1, ...],
                'L2_error': {
                    'field_name_1': [error_0, error_1, ...],
                    'field_name_2': [error_0, error_1, ...],
                    ...
                }
            }
            Onde 'field_name_1' poderia ser 'Ez', por exemplo.

        Retorna
        -------
        None
            A função exibe um gráfico e não retorna nenhum valor.
        """
        time_values = error_data['time']
        l2_error_dict = error_data['L2_error']

        # 1. Criar a figura e os eixos para o plot
        # O uso de subplots é uma boa prática para ter mais controle sobre a figura.
        _, ax = plt.subplots(figsize=(10, 6), dpi=90)

        # 2. Iterar sobre cada campo no dicionário de erros e plotar
        for field_name, errors in l2_error_dict.items():
            # Usamos semilogy para que o eixo y fique em escala logarítmica.
            # Isso é ideal para visualizar ordens de magnitude do erro.
            ax.semilogy(time_values, errors, label=f'Erro L2 - {field_name}', marker='o', markersize=3, linestyle='-')

        # 3. Configurar os detalhes do gráfico para clareza
        ax.set_title('Evolução do Erro $L^2$ ao Longo do Tempo', fontsize=16)
        ax.set_xlabel('Tempo (s)', fontsize=12)
        ax.set_ylabel('Erro $L^2$ (Escala Log)', fontsize=12)
        ax.legend()  # Adiciona a legenda (baseada nos 'labels' definidos no plot)
        ax.grid(True, which='both', linestyle='--', linewidth=0.5) # Adiciona uma grade para facilitar a leitura

        # 4. Ajustar o layout e exibir o gráfico
        plt.tight_layout()


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


class UniaxialPml:
    """
    Classe para simular e visualizar a propagação de ondas em 3D,
    incluindo a fonte de pulso gaussiano e métodos de plotagem.
    """
    def __init__(self, problem):
        self.problem = problem
        self.DIM = problem['DIM']
        self.Lx = problem['domain']['Lx']
        self.t_final = problem['dg']['t_final']        
        self.x0 = self.Lx - problem['pml']['L']


    def gaussian_pulse_3d(self, x, y, z, t, c=1.0):
        """
        Pulso gaussiano esfericamente simétrico deslocado e propagante no tempo em 3D.

        O pulso se origina em (x0, y0, z0).

        Ψ(x, y, z, t) = A * exp(- (sqrt((x-x₀)² + (y-y₀)² + (z-z₀)²) - c*t)² / (2σ²))

        Parâmetros:
        - x, y, z: coordenadas espaciais do ponto de avaliação
        - t: instante de tempo
        - x0, y0, z0: coordenadas do centro de origem do pulso
        - sigma: largura (espessura) da casca do pulso
        - amplitude: valor máximo do campo
        - c: velocidade de propagação

        Retorna:
        - Ψ(x, y, z, t): campo escalar em 3D
        """
        sigma_r = self.problem['source']['sigma_r']     # Desvio padrão espacial do pulso
        a0 = self.problem['source']['A0']               # Amplitude do pulso
        t0 = self.problem['source']['t0']               # Tempo de pico do pulso

        # Coordenadas do centro do pulso
        x0, y0, z0 = self.problem['source']['x0']
        
        # Calcula a distância radial 'r' do novo centro (x0, y0, z0) ao ponto (x, y, z)
        r = np.sqrt((x - x0)**2 + (y - y0)**2 + (z - z0)**2)
        
        # A fórmula do pulso é a mesma, mas agora baseada na distância 'r' deslocada
        return a0 * np.exp(-((r - c * (t - t0))**2) / (2 * sigma_r**2))


    def plot_L2_error(self, error_data) -> None:
        """
        Plota a evolução do erro L2 ao longo do tempo para um ou mais campos.

        A função gera um gráfico com o tempo no eixo x e o erro L2 no eixo y,
        utilizando uma escala logarítmica para melhor visualização de pequenas
        variações de erro.

        Parâmetros
        ----------
        error_data : dict
            Dicionário contendo os dados do erro. Espera-se a seguinte estrutura:
            {
                'time': [t_0, t_1, ...],
                'L2_error': {
                    'field_name_1': [error_0, error_1, ...],
                    'field_name_2': [error_0, error_1, ...],
                    ...
                }
            }
            Onde 'field_name_1' poderia ser 'Ez', por exemplo.

        Retorna
        -------
        None
            A função exibe um gráfico e não retorna nenhum valor.
        """
        time_values = error_data['time']
        l2_error_dict = error_data['L2_error']

        # 1. Criar a figura e os eixos para o plot
        # O uso de subplots é uma boa prática para ter mais controle sobre a figura.
        _, ax = plt.subplots(figsize=(10, 6), dpi=90)

        # 2. Iterar sobre cada campo no dicionário de erros e plotar
        for field_name, errors in l2_error_dict.items():
            # Usamos semilogy para que o eixo y fique em escala logarítmica.
            # Isso é ideal para visualizar ordens de magnitude do erro.
            ax.semilogy(time_values, errors, label=f'Erro L2 - {field_name}', marker='o', markersize=3, linestyle='-')

        # 3. Configurar os detalhes do gráfico para clareza
        ax.set_title('Evolução do Erro $L^2$ ao Longo do Tempo', fontsize=16)
        ax.set_xlabel('Tempo (s)', fontsize=12)
        ax.set_ylabel('Erro $L^2$ (Escala Log)', fontsize=12)
        ax.legend()  # Adiciona a legenda (baseada nos 'labels' definidos no plot)
        ax.grid(True, which='both', linestyle='--', linewidth=0.5) # Adiciona uma grade para facilitar a leitura

        # 4. Ajustar o layout e exibir o gráfico
        plt.tight_layout()

    
    def animate_gaussian_pulse_propagation(self, grid_size=150, slice_z=0.0, num_frames=150, output_filename=None):
        """
        Cria e exibe (ou salva) uma animação da propagação do pulso gaussiano.

        A animação mostra um corte transversal 2D do campo 3D no plano z especificado.

        Parâmetros
        ----------
        grid_size : int, opcional
            Número de pontos no grid de visualização (em cada direção).
        slice_z : float, opcional
            A coordenada z do plano de corte para a visualização.
        num_frames : int, opcional
            O número de quadros na animação.
        output_filename : str, opcional
            Se fornecido (ex: 'pulso.gif' ou 'pulso.mp4'), salva a animação 
            neste arquivo em vez de exibi-la na tela.
        """
        fig, ax = plt.subplots(figsize=(8, 7))

        # Define a grade espacial baseada no domínio do problema (self.Lx)
        x_lim = self.Lx
        x = np.linspace(-x_lim, x_lim, grid_size)
        y = np.linspace(-x_lim, x_lim, grid_size)
        X, Y = np.meshgrid(x, y)

        # Configura a normalização da barra de cores para manter a consistência
        norm = plt.Normalize(vmin=0, vmax=1.0)
        
        # Cria um plot inicial com a barra de cores
        quadro_inicial = np.zeros_like(X)
        cax = ax.pcolormesh(X, Y, quadro_inicial, cmap='viridis', norm=norm)
        fig.colorbar(cax, ax=ax, label="Amplitude do Campo")

        # Função interna que atualiza cada quadro da animação
        def update(frame):
            # Calcula o tempo atual baseado no número de quadros e no tempo final
            t = (frame / num_frames) * self.t_final
            
            # Calcula o valor do pulso gaussiano no corte 2D
            pulse_data = self.gaussian_pulse_3d(X, Y, slice_z, t)
            
            # Atualiza os dados do plot
            cax.set_array(pulse_data.ravel())
            
            ax.set_title(f'Propagação do Pulso 3D (Corte em z={slice_z:.2f}, t={t:.2f} s)')
            ax.set_xlabel('X')
            ax.set_ylabel('Y')
            ax.set_aspect('equal', 'box')
            
            return [cax]

        # Cria o objeto de animação
        ani = FuncAnimation(fig, update, frames=num_frames, blit=False, interval=40)

        # Decide se salva ou mostra a animação
        if output_filename:
            try:
                # Requer ffmpeg (para .mp4) ou imagemagick (para .gif)
                writer = 'ffmpeg' if output_filename.endswith('.mp4') else 'imagemagick'
                print(f"Salvando animação em '{output_filename}'... (Isso pode levar um momento)")
                ani.save(output_filename, writer=writer, fps=25)
                print(f"Animação salva com sucesso em '{output_filename}'.")
            except Exception as e:
                print(f"Erro ao salvar a animação: {e}")
                print("Verifique se você tem 'ffmpeg' ou 'imagemagick' instalado e acessível no seu PATH.")
            plt.close(fig) # Fecha a figura para não exibi-la
        else:
            plt.show()


    def plot_field(self, xi, yi, zi, field_data, title, slice_tol=0.05, cbar_label="Amplitude do Campo (Ez)"):
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

        # Adiciona uma condição para não plotar pontos onde x < 0.
        # Cria uma máscara booleana para todos os pontos onde x é maior ou igual a 0.
        maskPos_x = xi >= 0

        # Aplica a máscara a todas as coordenadas e aos dados do campo.
        # Apenas os dados que correspondem à condição (x >= 0) serão plotados.
        scatter = ax1.scatter(
            xi[maskPos_x], 
            yi[maskPos_x], 
            zi[maskPos_x], 
            c=field_data[maskPos_x], 
            cmap='viridis', 
            s=5, 
            alpha=0.6
        )


        # Usa 'field_data' para a cor dos pontos
        cbar = fig.colorbar(scatter, ax=ax1, shrink=0.7, aspect=15, label=cbar_label)
        ax1.set_title("Visualização 3D Completa")
        ax1.set_xlabel('X')
        ax1.set_ylabel('Y')
        ax1.set_zlabel('Z')
        ax1.view_init(elev=25, azim=-135)

        # --- Subplot 2: Corte no plano z=0 ---
        ax2 = fig.add_subplot(1, 2, 2)

        # Filtra os pontos que estão próximos ao plano z=0
        mask = np.abs(zi) < slice_tol
        x_slice, y_slice, field_slice = xi[mask], yi[mask], field_data[mask]

        # Verifica se há pontos suficientes para plotar
        if len(x_slice) > 3:
            # tricontourf é excelente para dados não estruturados
            contour = ax2.tricontourf(x_slice, y_slice, field_slice, levels=20, cmap='viridis')
            fig.colorbar(contour, ax=ax2, label=cbar_label)
            ax2.set_title(f"Corte da Solução no Plano z ≈ 0 (±{slice_tol})")
        else:
            ax2.text(0.5, 0.5, 'Nenhum dado encontrado para o corte.', ha='center', va='center')
            ax2.set_title("Corte da Solução no Plano z ≈ 0")

        # Coordenadas dos vértices do quadrado para fechar o loop
        x0 = self.x0
        x_square = [-x0, x0, x0, -x0, -x0]
        y_square = [x0, x0, -x0, -x0, x0]

        # Desenha o quadrado no gráfico 2D
        ax2.plot(x_square, y_square, color='red', linestyle='--', linewidth=2, label=f'Interface (Lado={self.x0})')

        ax2.set_xlabel('X')
        ax2.set_ylabel('Y')
        ax2.set_aspect('equal', adjustable='box')
        ax2.grid(True, linestyle='--', alpha=0.5)

        plt.tight_layout(rect=[0, 0, 1, 0.96])
        

    def plot_sigma_fields(self, x, y, z, sigma_x, sigma_y, sigma_z, title):
        """
        Plota a distribuição espacial 3D dos campos de condutividade da PML.

        Esta função cria três subplots 3D, cada um mostrando a malha de pontos
        (x, y, z) e usando a cor para representar a magnitude de sigma_x, sigma_y
        e sigma_z, respectivamente. Ideal para visualizar dados de malhas 2D/3D.

        Parâmetros
        ----------
        self : object
            A instância da classe que contém os parâmetros da simulação.
        x, y, z : np.ndarray
            Arrays de coordenadas dos nós da malha (podem ser multidimensionais).
        sigma_x, sigma_y, sigma_z : np.ndarray
            Arrays com os valores dos campos de condutividade nos nós da malha.
        title : str
            Título geral para a figura.
        """
        # Garante que os arrays de entrada sejam ' achatados' para a plotagem,
        # o que funciona para qualquer dimensão de entrada (1D, 2D, etc.).
        x_flat = x.flatten()
        y_flat = y.flatten()
        z_flat = z.flatten()
        sigma_x_flat = sigma_x.flatten()
        sigma_y_flat = sigma_y.flatten()
        sigma_z_flat = sigma_z.flatten()

        # Verificação de consistência do tamanho após o achatamento
        msg_template = "Incompatibilidade de tamanho após achatamento entre as coordenadas (tamanho {len_coord}) e o sigma (tamanho {len_sigma})"
        assert len(x_flat) == len(sigma_x_flat), msg_template.format(len_coord=len(x_flat), len_sigma=len(sigma_x_flat))
        assert len(y_flat) == len(sigma_y_flat), msg_template.format(len_coord=len(y_flat), len_sigma=len(sigma_y_flat))
        assert len(z_flat) == len(sigma_z_flat), msg_template.format(len_coord=len(z_flat), len_sigma=len(sigma_z_flat))

        fig = plt.figure(figsize=(24, 8))
        fig.suptitle(title, fontsize=20)

        # --- Subplot 1: Distribuição 3D de σx ---
        ax1 = fig.add_subplot(1, 3, 1, projection='3d')
        sc1 = ax1.scatter(x_flat, y_flat, z_flat, c=sigma_x_flat, cmap='gray_r', s=8, alpha=0.8,
                        vmin=0, vmax=max(sigma_x_flat.max(), sigma_y_flat.max(), sigma_z_flat.max()) if sigma_x_flat.size > 0 else 1)

        ax1.set_title(r"$\sigma_x(x)$", fontsize=16)
        ax1.set_xlabel("X")
        ax1.set_ylabel("Y")
        ax1.set_zlabel("Z")
        cbar1 = fig.colorbar(sc1, ax=ax1, shrink=0.6, aspect=20, pad=0.1)
        cbar1.set_label(r"$\sigma_x$")
        ax1.view_init(elev=25, azim=-135) # Ajusta o ângulo de visão

        # --- Subplot 2: Distribuição 3D de σy ---
        ax2 = fig.add_subplot(1, 3, 2, projection='3d')
        sc2 = ax2.scatter(x_flat, y_flat, z_flat, c=sigma_y_flat, cmap='gray_r', s=8, alpha=0.8,
                        vmin=0, vmax=max(sigma_x_flat.max(), sigma_y_flat.max(), sigma_z_flat.max()) if sigma_y_flat.size > 0 else 1)

        ax2.set_title(r"$\sigma_y(y)$", fontsize=16)
        ax2.set_xlabel("X")
        ax2.set_ylabel("Y")
        ax2.set_zlabel("Z")
        cbar2 = fig.colorbar(sc2, ax=ax2, shrink=0.6, aspect=20, pad=0.1)
        cbar2.set_label(r"$\sigma_y$")
        ax2.view_init(elev=25, azim=-135)

        # --- Subplot 3: Distribuição 3D de σz ---
        ax3 = fig.add_subplot(1, 3, 3, projection='3d')
        sc3 = ax3.scatter(x_flat, y_flat, z_flat, c=sigma_z_flat, cmap='gray_r', s=8, alpha=0.8,
                        vmin=0, vmax=max(sigma_x_flat.max(), sigma_y_flat.max(), sigma_z_flat.max()) if sigma_z_flat.size > 0 else 1)
        
        ax3.set_title(r"$\sigma_z(z)$", fontsize=16)
        ax3.set_xlabel("X")
        ax3.set_ylabel("Y")
        ax3.set_zlabel("Z")
        cbar3 = fig.colorbar(sc3, ax=ax3, shrink=0.6, aspect=20, pad=0.1)
        cbar3.set_label(r"$\sigma_z$")
        ax3.view_init(elev=25, azim=-135)

        plt.tight_layout(rect=[0, 0, 1, 0.93])


class PlaneWaveExcitation:
    """
    Classe para simular e visualizar a propagação de ondas em 3D,
    incluindo a fonte de pulso gaussiano e métodos de plotagem.
    """
    def __init__(self, problem):
        self.problem = problem
        self.DIM = problem['DIM']
        self.Lx = problem['domain']['Lx']
        self.xa = problem['domain']['xa']
        self.x0 = problem['domain']['x0']
        self.t_final = problem['dg']['t_final']

        self.wavelength = problem['source']['wavelength']   # Comprimento de onda
        self.kx = 2 * np.pi / self.wavelength               # Número de onda na direção x
        self.ky = 0.0                                       # Número de onda na direção y (0 para onda plana)
        self.kz = 0.0                                       # Número de onda na direção z (0 para onda plana)
        self.omega = self.kx * problem['source']['c']       # Frequência angular


    def incident_Ey(self, x, y, z, t):
        """
        Calcula o campo elétrico Ey de uma onda plana incidente.

        A onda é definida como:
        Ey_inc(x, y, z, t) = A * cos(ω * t - kx * x - ky * y - kz * z)


        Parâmetros
        ----------
        A : float
            Amplitude do campo elétrico em V/m.
        w : float
            Frequência angular da onda (ω = 2πf) em rad/s. 
        r : tuple
            Coordenadas espaciais (x, y, z) do ponto onde o campo é avaliado.
        t : float
            Tempo no qual o campo é avaliado.

        Retorna
        -------
        float
            Valor do campo elétrico Ey no ponto e tempo especificados.
        """
        E0 = self.problem['source']['E0']
        return E0 * np.cos(self.omega * t - self.kx * x - self.ky * y - self.kz * z)


    def incident_Hz(self, x, y, z, t):
        """
        Calcula o campo magnético Hz de uma onda plana incidente.

        A onda é definida como:
        Hz_inc(x, y, z, t) = (1/eta) * cos(ω * t - kx * x - ky * y - kz * z)

        Parâmetros
        ----------
        r : tuple
            Coordenadas espaciais (x, y, z) do ponto onde o campo é avaliado.
        t : float
            Tempo no qual o campo é avaliado.

        Retorna
        -------
        float
            Valor do campo magnético Hz no ponto e tempo especificados.
        """
        eta = 1 # Impedância característica do meio (normalizada)
        return (1 / eta) * np.cos(self.omega * t - self.kx * x - self.ky * y - self.kz * z)


    def get_incident_fields(self, x, y, z, t):
        """
        Calcula todos os 6 componentes de campo de uma onda plana incidente.
        Para uma onda se propagando em +x, os campos não nulos são Ey e Hz.
        """
        Ex_inc = np.zeros_like(x)
        Ey_inc = self.incident_Ey(x, y, z, t)
        Ez_inc = np.zeros_like(x)        
        Hx_inc = np.zeros_like(x)
        Hy_inc = np.zeros_like(x)
        Hz_inc = self.incident_Hz(x, y, z, t)

        return Hx_inc, Hy_inc, Hz_inc, Ex_inc, Ey_inc, Ez_inc
    

    # def plot_incident_fields(self, sp, t):
    #     """
    #     Gera uma visualização 3D dos campos incidentes Ey e Hz AVALIADOS
    #     EXCLUSIVAMENTE nos nós da interface TF/SF.

    #     Este método de depuração mostra os valores exatos da fonte nos locais
    #     onde ela é injetada na simulação.

    #     Parâmetros
    #     ----------
    #     sp : Maxwell3D
    #         A instância da classe de discretização que contém os mapas da interface.
    #     t : float
    #         O instante de tempo "t" no qual o snapshot dos campos será gerado.
    #     """
    #     print(f"🔎 Gerando snapshot 3D dos campos incidentes na interface em t = {t:.2f}s...")

    #     # 1. Obter as coordenadas dos nós da interface a partir do objeto 'sp'
    #     # Esta lógica é a mesma do método plot_interface_nodes
    #     sf_id = sp.sf_group_id
    #     map_name = f'mapI_G{sf_id}'

    #     if not hasattr(sp, map_name):
    #         print(f"⚠️ Aviso: Mapa de interface '{map_name}' não encontrado. A plotagem foi cancelada.")
    #         return

    #     node_map = getattr(sp, map_name)
    #     vmap_indices = sp.vmapM[node_map]
    #     x_int = sp.x.ravel('F')[vmap_indices]
    #     y_int = sp.y.ravel('F')[vmap_indices]
    #     z_int = sp.z.ravel('F')[vmap_indices]

    #     # 2. Calcular os campos incidentes APENAS nesses pontos da interface
    #     Hx, Hy, Hz, Ex, Ey, Ez = self.get_incident_fields(x_int, y_int, z_int, t)

    #     # 3. Criar a figura com dois subplots 3D
    #     fig, (ax1, ax2) = plt.subplots(
    #         1, 2,
    #         figsize=(18, 8),
    #         subplot_kw={'projection': '3d'}
    #     )
    #     fig.suptitle(f'Campos Incidentes na Interface TF/SF em t = {t:.2f} s', fontsize=18)

    #     # --- Subplot para o Campo Elétrico Ey na interface ---
    #     scatter1 = ax1.scatter(x_int, y_int, z_int, c=Ey, cmap='viridis', s=25, depthshade=True)
    #     ax1.set_title(r"Campo Elétrico $E_y$ na Interface")
    #     fig.colorbar(scatter1, ax=ax1, shrink=0.6, aspect=20, label=r'Amplitude $E_y$ (V/m)')

    #     # --- Subplot para o Campo Magnético Hz na interface ---
    #     scatter2 = ax2.scatter(x_int, y_int, z_int, c=Hz, cmap='plasma', s=25, depthshade=True)
    #     ax2.set_title(r"Campo Magnético $H_z$ na Interface")
    #     fig.colorbar(scatter2, ax=ax2, shrink=0.6, aspect=20, label=r'Amplitude $H_z$ (A/m)')

    #     # 4. Configurações comuns para ambos os gráficos
    #     for ax in [ax1, ax2]:
    #         # Plota todos os nós da malha ao fundo para dar contexto
    #         ax.scatter(sp.x, sp.y, sp.z, c='gray', s=1, alpha=0.05)
    #         ax.set_xlabel("Eixo X")
    #         ax.set_ylabel("Eixo Y")
    #         ax.set_zlabel("Eixo Z")
    #         ax.view_init(elev=30, azim=-60)

    #     plt.tight_layout(rect=[0, 0, 1, 0.95])

    def plot_incident_fields(self, sp, group_id_A, group_id_B, t):
        """
        Gera uma visualização 3D dos campos incidentes Ey e Hz AVALIADOS
        EXCLUSIVAMENTE nos nós de uma interface específica.

        Este método de depuração mostra os valores exatos da fonte nos locais
        onde ela é injetada na simulação, para a interface entre os grupos A e B.

        Parâmetros
        ----------
        sp : Maxwell3D
            A instância da classe de discretização que contém os mapas da interface.
        group_id_A, group_id_B : int
            Os IDs dos dois grupos que formam a interface a ser visualizada.
        t : float
            O instante de tempo "t" no qual o snapshot dos campos será gerado.
        """
        # 1. Obter nomes descritivos dos grupos a partir do objeto de malha
        name_A = sp.mesh.group_names.get(group_id_A, f"ID {group_id_A}")
        name_B = sp.mesh.group_names.get(group_id_B, f"ID {group_id_B}")
        print(f"🔎 Gerando snapshot dos campos incidentes na interface entre '{name_A}' e '{name_B}' em t = {t:.2f}s...")

        # 2. Obter os IDs dos nós da interface a partir da perspectiva do grupo A
        try:
            # vmapIM_G{id_A} contém os IDs globais dos nós no lado 'Minus' (grupo A) da interface
            node_ids = getattr(sp, f'vmapIM_G{group_id_A}')
        except AttributeError:
            print(f"⚠️ Aviso: Mapa de interface 'vmapIM_G{group_id_A}' não encontrado.")
            print("-> Certifique-se de que 'buildGroupInterfaceMaps' foi executado e a interface foi descoberta.")
            return

        if node_ids.size == 0:
            print(f"⚠️ Nenhum nó de interface encontrado para o grupo '{name_A}' na fronteira com '{name_B}'.")
            return

        # 3. Obter as coordenadas (x,y,z) desses nós específicos
        x_int = sp.x.ravel('F')[node_ids]
        y_int = sp.y.ravel('F')[node_ids]
        z_int = sp.z.ravel('F')[node_ids]

        # 4. Calcular os campos incidentes APENAS nesses pontos da interface
        Hx, Hy, Hz, Ex, Ey, Ez = self.get_incident_fields(x_int, y_int, z_int, t)

        # 5. Criar a figura com dois subplots 3D
        fig, (ax1, ax2) = plt.subplots(
            1, 2,
            figsize=(18, 8),
            subplot_kw={'projection': '3d'}
        )
        fig.suptitle(f"Campos Incidentes na Interface ('{name_A}' / '{name_B}') em t = {t:.2f} s", fontsize=18)

        # --- Subplot para o Campo Elétrico Ey na interface ---
        scatter1 = ax1.scatter(x_int, y_int, z_int, c=Ey, cmap='viridis', s=25, depthshade=True)
        ax1.set_title(r"Campo Elétrico $E_y$ na Interface")
        fig.colorbar(scatter1, ax=ax1, shrink=0.6, aspect=20, label=r'Amplitude $E_y$ (V/m)')

        # --- Subplot para o Campo Magnético Hz na interface ---
        scatter2 = ax2.scatter(x_int, y_int, z_int, c=Hz, cmap='plasma', s=25, depthshade=True)
        ax2.set_title(r"Campo Magnético $H_z$ na Interface")
        fig.colorbar(scatter2, ax=ax2, shrink=0.6, aspect=20, label=r'Amplitude $H_z$ (A/m)')

        # 6. Configurações comuns para ambos os gráficos
        for ax in [ax1, ax2]:
            ax.scatter(sp.x, sp.y, sp.z, c='gray', s=1, alpha=0.05) # Nós da malha ao fundo
            ax.set_xlabel("Eixo X")
            ax.set_ylabel("Eixo Y")
            ax.set_zlabel("Eixo Z")
            ax.view_init(elev=30, azim=-60)

        plt.tight_layout(rect=[0, 0, 1, 0.95])

    def plot_L2_error(self, error_data) -> None:
        """
        Plota a evolução do erro L2 ao longo do tempo para um ou mais campos.

        A função gera um gráfico com o tempo no eixo x e o erro L2 no eixo y,
        utilizando uma escala logarítmica para melhor visualização de pequenas
        variações de erro.

        Parâmetros
        ----------
        error_data : dict
            Dicionário contendo os dados do erro. Espera-se a seguinte estrutura:
            {
                'time': [t_0, t_1, ...],
                'L2_error': {
                    'field_name_1': [error_0, error_1, ...],
                    'field_name_2': [error_0, error_1, ...],
                    ...
                }
            }
            Onde 'field_name_1' poderia ser 'Ez', por exemplo.

        Retorna
        -------
        None
            A função exibe um gráfico e não retorna nenhum valor.
        """
        time_values = error_data['time']
        l2_error_dict = error_data['L2_error']

        # 1. Criar a figura e os eixos para o plot
        # O uso de subplots é uma boa prática para ter mais controle sobre a figura.
        _, ax = plt.subplots(figsize=(10, 6), dpi=90)

        # 2. Iterar sobre cada campo no dicionário de erros e plotar
        for field_name, errors in l2_error_dict.items():
            # Usamos semilogy para que o eixo y fique em escala logarítmica.
            # Isso é ideal para visualizar ordens de magnitude do erro.
            ax.semilogy(time_values, errors, label=f'Erro L2 - {field_name}', marker='o', markersize=3, linestyle='-')

        # 3. Configurar os detalhes do gráfico para clareza
        ax.set_title('Evolução do Erro $L^2$ ao Longo do Tempo', fontsize=16)
        ax.set_xlabel('Tempo (s)', fontsize=12)
        ax.set_ylabel('Erro $L^2$ (Escala Log)', fontsize=12)
        ax.legend()  # Adiciona a legenda (baseada nos 'labels' definidos no plot)
        ax.grid(True, which='both', linestyle='--', linewidth=0.5) # Adiciona uma grade para facilitar a leitura

        # 4. Ajustar o layout e exibir o gráfico
        plt.tight_layout()


    def plot_field(self, xi, yi, zi, field_data, title, slice_tol=0.05, cbar_label="Amplitude do Campo (Ez)"):
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

        # Adiciona uma condição para não plotar pontos onde x < 0.
        # Cria uma máscara booleana para todos os pontos onde x é maior ou igual a 0.
        maskPos_x = xi >= 0

        # Aplica a máscara a todas as coordenadas e aos dados do campo.
        # Apenas os dados que correspondem à condição (x >= 0) serão plotados.
        scatter = ax1.scatter(
            xi[maskPos_x], 
            yi[maskPos_x], 
            zi[maskPos_x], 
            c=field_data[maskPos_x], 
            cmap='viridis', 
            s=5, 
            alpha=0.6
        )


        # Usa 'field_data' para a cor dos pontos
        cbar = fig.colorbar(scatter, ax=ax1, shrink=0.7, aspect=15, label=cbar_label)
        ax1.set_title("Visualização 3D Completa")
        ax1.set_xlabel('X')
        ax1.set_ylabel('Y')
        ax1.set_zlabel('Z')
        ax1.view_init(elev=25, azim=-135)

        # --- Subplot 2: Corte no plano z=0 ---
        ax2 = fig.add_subplot(1, 2, 2)

        # Filtra os pontos que estão próximos ao plano z=0
        mask = np.abs(zi) < slice_tol
        x_slice, y_slice, field_slice = xi[mask], yi[mask], field_data[mask]

        # Verifica se há pontos suficientes para plotar
        if len(x_slice) > 3:
            # tricontourf é excelente para dados não estruturados
            contour = ax2.tricontourf(x_slice, y_slice, field_slice, levels=20, cmap='viridis')
            fig.colorbar(contour, ax=ax2, label=cbar_label)
            ax2.set_title(f"Corte da Solução no Plano z ≈ 0 (±{slice_tol})")
        else:
            ax2.text(0.5, 0.5, 'Nenhum dado encontrado para o corte.', ha='center', va='center')
            ax2.set_title("Corte da Solução no Plano z ≈ 0")

        # Coordenadas dos vértices do quadrado para fechar o loop
        x0 = self.x0
        x_square = [-x0, x0, x0, -x0, -x0]
        y_square = [x0, x0, -x0, -x0, x0]

        # Desenha o quadrado no gráfico 2D
        ax2.plot(x_square, y_square, color='red', linestyle='--', linewidth=2, label=f'Interface (Lado={self.x0})')

        ax2.set_xlabel('X')
        ax2.set_ylabel('Y')
        ax2.set_aspect('equal', adjustable='box')
        ax2.grid(True, linestyle='--', alpha=0.5)

        plt.tight_layout(rect=[0, 0, 1, 0.96])
        

    def plot_sigma_fields(self, x, y, z, sigma_x, sigma_y, sigma_z, title):
        """
        Plota a distribuição espacial 3D dos campos de condutividade da PML.

        Esta função cria três subplots 3D, cada um mostrando a malha de pontos
        (x, y, z) e usando a cor para representar a magnitude de sigma_x, sigma_y
        e sigma_z, respectivamente. Ideal para visualizar dados de malhas 2D/3D.

        Parâmetros
        ----------
        self : object
            A instância da classe que contém os parâmetros da simulação.
        x, y, z : np.ndarray
            Arrays de coordenadas dos nós da malha (podem ser multidimensionais).
        sigma_x, sigma_y, sigma_z : np.ndarray
            Arrays com os valores dos campos de condutividade nos nós da malha.
        title : str
            Título geral para a figura.
        """
        # Garante que os arrays de entrada sejam ' achatados' para a plotagem,
        # o que funciona para qualquer dimensão de entrada (1D, 2D, etc.).
        x_flat = x.flatten()
        y_flat = y.flatten()
        z_flat = z.flatten()
        sigma_x_flat = sigma_x.flatten()
        sigma_y_flat = sigma_y.flatten()
        sigma_z_flat = sigma_z.flatten()

        # Verificação de consistência do tamanho após o achatamento
        msg_template = "Incompatibilidade de tamanho após achatamento entre as coordenadas (tamanho {len_coord}) e o sigma (tamanho {len_sigma})"
        assert len(x_flat) == len(sigma_x_flat), msg_template.format(len_coord=len(x_flat), len_sigma=len(sigma_x_flat))
        assert len(y_flat) == len(sigma_y_flat), msg_template.format(len_coord=len(y_flat), len_sigma=len(sigma_y_flat))
        assert len(z_flat) == len(sigma_z_flat), msg_template.format(len_coord=len(z_flat), len_sigma=len(sigma_z_flat))

        fig = plt.figure(figsize=(24, 8))
        fig.suptitle(title, fontsize=20)

        # --- Subplot 1: Distribuição 3D de σx ---
        ax1 = fig.add_subplot(1, 3, 1, projection='3d')
        sc1 = ax1.scatter(x_flat, y_flat, z_flat, c=sigma_x_flat, cmap='gray_r', s=8, alpha=0.8,
                        vmin=0, vmax=max(sigma_x_flat.max(), sigma_y_flat.max(), sigma_z_flat.max()) if sigma_x_flat.size > 0 else 1)

        ax1.set_title(r"$\sigma_x(x)$", fontsize=16)
        ax1.set_xlabel("X")
        ax1.set_ylabel("Y")
        ax1.set_zlabel("Z")
        cbar1 = fig.colorbar(sc1, ax=ax1, shrink=0.6, aspect=20, pad=0.1)
        cbar1.set_label(r"$\sigma_x$")
        ax1.view_init(elev=25, azim=-135) # Ajusta o ângulo de visão

        # --- Subplot 2: Distribuição 3D de σy ---
        ax2 = fig.add_subplot(1, 3, 2, projection='3d')
        sc2 = ax2.scatter(x_flat, y_flat, z_flat, c=sigma_y_flat, cmap='gray_r', s=8, alpha=0.8,
                        vmin=0, vmax=max(sigma_x_flat.max(), sigma_y_flat.max(), sigma_z_flat.max()) if sigma_y_flat.size > 0 else 1)

        ax2.set_title(r"$\sigma_y(y)$", fontsize=16)
        ax2.set_xlabel("X")
        ax2.set_ylabel("Y")
        ax2.set_zlabel("Z")
        cbar2 = fig.colorbar(sc2, ax=ax2, shrink=0.6, aspect=20, pad=0.1)
        cbar2.set_label(r"$\sigma_y$")
        ax2.view_init(elev=25, azim=-135)

        # --- Subplot 3: Distribuição 3D de σz ---
        ax3 = fig.add_subplot(1, 3, 3, projection='3d')
        sc3 = ax3.scatter(x_flat, y_flat, z_flat, c=sigma_z_flat, cmap='gray_r', s=8, alpha=0.8,
                        vmin=0, vmax=max(sigma_x_flat.max(), sigma_y_flat.max(), sigma_z_flat.max()) if sigma_z_flat.size > 0 else 1)
        
        ax3.set_title(r"$\sigma_z(z)$", fontsize=16)
        ax3.set_xlabel("X")
        ax3.set_ylabel("Y")
        ax3.set_zlabel("Z")
        cbar3 = fig.colorbar(sc3, ax=ax3, shrink=0.6, aspect=20, pad=0.1)
        cbar3.set_label(r"$\sigma_z$")
        ax3.view_init(elev=25, azim=-135)

        plt.tight_layout(rect=[0, 0, 1, 0.93])

