# FM235 - Dinâmica de Missões Espaciais Modernas
# Author: Mairon de Souza Wolniewicz
# Date: 2025-Sep-10

# 6. Análise de estabilidade linear dos equilíbrios do PR3C.

# Parte I.
# Calcule os autovalores da matriz jacobiana associada aos equilíbrios L1 a L5, para 1e−6 < µ < 0.5.
# Apresente graficamente estes resultados.

# Parte II.
# Para os sistemas (i) Terra-Lua e (ii) Sol-Terra, apresente os valores numéricos dos 6 auto-valores de cada equilíbrio.
# Classifique cada ponto de equilíbrio e comente os aspectos de estabilidade linear.
# Em particular, compare L3 com L1 e L2.

import utils
import numpy as np
import matplotlib.pyplot as plt

def plot_equilibirium_eigenvalues(eigenvalues_dict, title=""):
    # Define colors for the plot
    color_list = ["tab:blue","tab:orange","tab:green","tab:red","tab:purple","tab:brown"]
    # Create plot
    fig, axes = plt.subplots(2,3)
    ax = axes.flatten()
    fig.delaxes(ax[-1])
    ax = ax[:-1]

    i = 0
    for L_case in eigenvalues_dict.keys():
        for j in range(6):
            # Compute labels
            m = j // 2 + 1 # increases every 2 steps
            n = j % 2 + 1  # alternates between 1 and 2
            ax[i].plot(eigenvalues_dict[L_case][j].real,     eigenvalues_dict[L_case][j].imag,     '.', color=color_list[j])
            ax[i].plot(eigenvalues_dict[L_case][j][0].real,  eigenvalues_dict[L_case][j][0].imag,  'x', color=color_list[j], label=f"$\lambda_{m},_{n}$")
            ax[i].plot(eigenvalues_dict[L_case][j][-1].real, eigenvalues_dict[L_case][j][-1].imag, 'x', color=color_list[j])

        ax[i].set_title(f"{title} {L_case}")
        ax[i].axvline(0, color='k', linestyle='-', linewidth=0.8)
        ax[i].axhline(0, color='k', linestyle='-', linewidth=0.8)
        ax[i].set_xlabel('Re($\lambda$)')
        ax[i].set_ylabel('Im($\lambda$)')
        ax[i].grid()
        ax[i].legend()
        i+=1

# For a list of µ, compute the eigenvalues
def compute_eigenvalues(mu_list, convention="Barcelona"):
    n_iters = len(mu_list)
    eigenvalues_dict = {
        "L1" : np.zeros((6,n_iters),dtype=np.complex128),
        "L2" : np.zeros((6,n_iters),dtype=np.complex128),
        "L3" : np.zeros((6,n_iters),dtype=np.complex128),
        "L4" : np.zeros((6,n_iters),dtype=np.complex128),
        "L5" : np.zeros((6,n_iters),dtype=np.complex128),
    }

    for i, mu in enumerate(mu_list):
        # Compute equilibrium points coordinates
        equilibrium_points = utils.equilibrium_points(mu, convention=convention)
        for L_case in eigenvalues_dict.keys():
            # For each equilibrium point, compute the eigenvalues
            x = equilibrium_points[L_case][0]
            y = equilibrium_points[L_case][1]
            z = 0
            Dxf = utils.linearization_matrix(x, y, z, mu, convention=convention)
            eigenvalues = np.linalg.eigvals(Dxf)

            # Store the eigenvalues for plotting
            for j in range(6):
                eigenvalues_dict[L_case][j,i] = eigenvalues[j]
    return eigenvalues_dict

# Define primary position convention
convention = "Barcelona"

# Mass of celestial bodies [kg]
body_mass = {
    "Sun"   : 1.988500e+30,
    "Earth" : 5.972190e+24,
    "Moon"  : 7.349000e+22,
}

# Compute the mass parameter µ for each system
system_mu = {
    "Earth-Moon" : utils.mu_from_masses(body_mass["Earth"], body_mass["Moon"]),
    "Sun-Earth"  : utils.mu_from_masses(body_mass["Sun"],   body_mass["Earth"]),
}

# Variable µ list
mu_list = np.linspace(1e-6, 0.5, 100)
eigenvalues_dict = compute_eigenvalues(mu_list, convention=convention)
plot_equilibirium_eigenvalues(eigenvalues_dict, "1e-6 < µ < 0.5")

# Earth-Moon
equilibrium_points = utils.equilibrium_points(system_mu["Earth-Moon"], convention=convention)
eigenvalues_dict   = compute_eigenvalues([system_mu["Earth-Moon"]])
plot_equilibirium_eigenvalues(eigenvalues_dict, "Earth-Moon")
print("Eigenvalues for Earth-Moon System")
for L_case in eigenvalues_dict.keys():
    print(L_case)
    print(eigenvalues_dict[L_case])

# Sun-Earth
equilibrium_points = utils.equilibrium_points(system_mu["Sun-Earth"], convention=convention)
eigenvalues_dict   = compute_eigenvalues([system_mu["Sun-Earth"]])
plot_equilibirium_eigenvalues(eigenvalues_dict, "Sun-Earth")
print("Eigenvalues for Sun-Earth System")
for L_case in eigenvalues_dict.keys():
    print(L_case)
    print(eigenvalues_dict[L_case])

plt.show()
