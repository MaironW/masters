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

def linearization_matrix(x, y, z, mu, convention="Barcelona"):
    # Build the equilibrium matrix Dxf
    Dxf = np.zeros([6,6])

    Dxf[0,3] =  1 # df1dx4
    Dxf[1,4] =  1 # df2dx5
    Dxf[2,5] =  1 # df3dx6
    Dxf[3,4] =  2 # df4dx5
    Dxf[4,3] = -2 # df5dx4

    # Barcelona parameters
    A = x - mu
    B = x - mu + 1
    r1 = ((x - mu)**2 + y**2 + z**2)**0.5
    r2 = ((x - mu + 1)**2 + y**2 + z**2)**0.5

    # Compute the Hessian
    Uxx = 1 - (1-mu)/r1**3 - mu/r2**3 + 3*(1-mu)*A**2/r1**5 + 3*mu*B**2/r2**5
    Uyy = 1 - (1-mu)/r1**3 - mu/r2**3 + 3*(1-mu)*y**2/r1**5 + 3*mu*y**2/r2**5
    Uzz =   - (1-mu)/r1**3 - mu/r2**3 + 3*(1-mu)*z**2/r1**5 + 3*mu*z**2/r2**5
    Uxy = 3*(1-mu)*A*y/r1**5 + 3*mu*B*y/r2**5
    Uxz = 3*(1-mu)*A*z/r1**5 + 3*mu*B*z/r2**5
    Uyz = 3*(1-mu)*y*z/r1**5 + 3*mu*y*z/r2**5

    Dxf[3:6, 0:3] = np.array([[Uxx, Uxy, Uxz],
                              [Uxy, Uyy, Uyz],
                              [Uxz, Uyz, Uzz]])
    return Dxf

n_iters = 100
mu_list = np.linspace(1e-6, 0.5, n_iters)

eigenvalues_dict = {
    "L1" : np.zeros((6,n_iters),dtype=np.complex128),
    "L2" : np.zeros((6,n_iters),dtype=np.complex128),
    "L3" : np.zeros((6,n_iters),dtype=np.complex128),
    "L4" : np.zeros((6,n_iters),dtype=np.complex128),
}

for i in range(n_iters):
    # Compute equilibrium points coordinates
    equilibrium_points = utils.equilibrium_points(mu_list[i], convention="Barcelona")
    for L_case in eigenvalues_dict.keys():

        # For each equilibrium point, compute the eigenvalues
        x = equilibrium_points[L_case][0]
        y = equilibrium_points[L_case][1]
        z = 0
        Dxf = linearization_matrix(x, y, z, mu_list[i], convention="Barcelona")
        eigenvalues = np.linalg.eigvals(Dxf)

        # Store the eigenvalues for plotting
        for j in range(6):
            eigenvalues_dict[L_case][j,i] = eigenvalues[j]

# Define colors for the plot
color_list = ["tab:blue","tab:orange","tab:green","tab:red","tab:purple","tab:brown"]

# Create plot
fig, axes = plt.subplots(1,4)
ax = axes.flatten()
i = 0
for L_case in eigenvalues_dict.keys():
    for j in range(6):
        ax[i].plot(eigenvalues_dict[L_case][j].real,     eigenvalues_dict[L_case][j].imag,     '.', color=color_list[j])
        ax[i].plot(eigenvalues_dict[L_case][j][0].real,  eigenvalues_dict[L_case][j][0].imag,  'x', color=color_list[j], label=f"$\lambda_{j+1}$")
        ax[i].plot(eigenvalues_dict[L_case][j][-1].real, eigenvalues_dict[L_case][j][-1].imag, 'o', color=color_list[j])

    ax[i].set_title(L_case)
    ax[i].axvline(0, color='gray', linestyle='--', linewidth=0.8)
    ax[i].axhline(0, color='gray', linestyle='--', linewidth=0.8)
    ax[i].set_xlabel('Real Part')
    ax[i].set_ylabel('Imaginary Part')
    ax[i].grid()
    ax[i].legend()
    i+=1
plt.show()