# FM235 - Dinâmica de Missões Espaciais Modernas
# Author: Mairon de Souza Wolniewicz
# Date: 2025-Nov-02
# Messy code, but computes the Lyapunov orbits around L1 of the Sun-Jupiter system

import utils
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from scipy.integrate import solve_ivp

###############
# DEFINITIONS #
###############

# Mass of celestial bodies [kg]
body_mass = {
    "Sun"     : 1.988500e30,
    "Earth"   : 5.972190e24,
    "Jupiter" : 1.898130e27,
    "Moon"    : 7.349000e22,
}

# Compute the mass parameter µ for the Sun-Jupiter system
mu = utils.mu_from_masses(body_mass["Sun"], body_mass["Jupiter"])

# Compute equilibrium points
equilibrium_points = utils.equilibrium_points(mu)

#############
# FUNCTIONS #
#############

# Generate the variational equations for the extended state
def variational_equations(t, state_M, mu):
    # Extract the state vector and the matrix M
    state  = state_M[:6]
    M_flat = state_M[6:]
    M      = M_flat.reshape((6,6))
    # Dynamic equations
    dstate = utils.state_equations(t, state, mu)
    # Jacobian matrix A(t)
    A = utils.linearization_matrix(state[0], state[1], state[2], mu)
    # Variational equations: dM/dt = A*M
    dM = A @ M
    dM_flat = dM.flatten()
    # Return the concatenated vector (6 + 36 variables)
    output = np.empty_like(state_M)
    output[:6] = dstate
    output[6:] = dM_flat
    return output

# Integrate the state transition matrix up to T/2, when the solution crosses y=0
def integrate_variational(mu, init_state, t_end, terminate=True):
    def y_event(t, y, mu):
        return y[1]
    y_event.terminal = terminate
    y_event.direction = -1

    # M =  dx/dx  dx/dy  dx/dz  dx/dvx   dx/dvy   dx/dvz
    #      dy/dx  dy/dy  dy/dz  dy/dvx  [dy/dvy]  dy/dvz
    #      dz/dx  dz/dy  dz/dz  dz/dvx   dz/dvy   dz/dvz
    #     dvx/dx dvx/dy dvx/dz dvx/dvx [dvx/dvy] dvx/dvz
    #     dvy/dx dvy/dy dvy/dz dvy/dvx  dvy/dvy  dvy/dvz
    #     dvz/dx dvz/dy dvz/dz dvz/dvx  dvz/dvy  dvz/dvz
    M0 = np.eye(6).ravel()
    y0 = np.concatenate((init_state, M0))
    t_span = (0, t_end)
    solution = solve_ivp(variational_equations, t_span, y0, args=(mu,), rtol=1e-13, atol=1e-15, events=y_event)
    return solution

# Single shoot algorithm to find the planar orbits
def differential_correction(mu, init_state, t_end, tol=1e-8, n_iter=10):
    # Make a copy of the initial state to not alter it
    state = np.copy(init_state)
    T_half = t_end
    for i in range(n_iter):
        # Integrate the augmented state
        solution = integrate_variational(mu, state, T_half, terminate=True)
        # The solution may be valid if the solution crosses the y=0 axis
        if solution.t_events[0].size > 0:
            t_end = solution.t_events[0][0]
            state_half = solution.y[:6,-1]
            deriv_half = utils.state_equations(t_end, state_half, mu)
            y_half  = state_half[1]
            vx_half = state_half[3]
            vy_half = state_half[4]
            ax_half = deriv_half[3]

            # Correction condition
            F = np.array([y_half, vx_half])
            if np.linalg.norm(F) < tol:
                break

            # Get useful parameters from the state transition matrix at half period
            M = solution.y[6:,-1].reshape((6,6))
            dy_dvy  = M[1,4]
            dvx_dvy = M[3,4]

            # Build the minimal state transition matrix
            M_sub = np.array([[ dy_dvy, vy_half],
                              [dvx_dvy, ax_half]])

            # Correction of [vy0 only]
            delta = np.linalg.solve(M_sub, -F)
            state[4] += delta[0]
    return solution

# Get solutions
def compute_solutions(x_ini, x_end, n=1, tol=1e-8):
    # Lists to store solutions
    init_state_list = []
    tf = 1.55
    vy = 0.01
    for x0 in np.linspace(x_ini, x_end, n):
        # Brute force first guesses for vy0
        vy0_new = None
        for vy0 in np.linspace(vy, 0.6, 25):
            # vy0_new will be None until a valid solution is found.
            # In this case, keep using the brute force guesses.
            if vy0_new is not None:
                vy0 = vy0_new
            init_state = np.array([x0, 0, 0, 0, vy0, 0])
            solution = differential_correction(mu, init_state, tf, tol=tol)

            if solution.t_events[0].size > 0:
                T_half = solution.t_events[0][0]
                x_ini  = solution.y[0,0]
                x_half = solution.y[0,-1]
                # If xf > x0, save updated init state
                if x_half > x_ini:
                    x,y,z,vx,vy,vz = solution.y[:6,0]
                    print(x,vy,T_half)
                    vy0_new = vy
                    init_state_list.append([x,y,z,vx,vy,vz,T_half])
                    tf = T_half*2
                    break # Stop iterating the velocity if solution found
        print(vy0_new)
        if vy0_new is None:
            break
    return init_state_list

# Load solutions from file in format [x,y,z,vx,vy,vz,T_half]
def load_solutions(file, n=None):
    init_state_list = np.loadtxt(file, delimiter=",", max_rows=n)
    return init_state_list

# Integrate a full orbit given initial states
def integrate_solutions(init_state_list):
    full_state_list = []
    for init_state in init_state_list:
        T_half = init_state[6]
        solution = integrate_variational(mu, init_state[:6], 2*T_half, terminate=False)
        full_state = solution.y
        full_state_list.append(full_state)
    return full_state_list

def stability_indices(M, tol=1e-6):
    M = np.asarray(M)
    eigvals = np.linalg.eigvals(M)
    s1 = 0.5*(eigvals[0] + 1/eigvals[0])
    s2 = 0.5*(eigvals[2] + 1/eigvals[2])
    s3 = 0.5*(eigvals[4] + 1/eigvals[4])
    eigvals[1] = 1/eigvals[0]
    eigvals[3] = 1/eigvals[2]
    eigvals[5] = 1/eigvals[4]
    return s1, s2, s3, eigvals

# Generic plot function to show results
def plot(x, y_sets, x_label=None, y_label=None, vlines=None, hlines=None, title=None, markers=None, type_list=None, legend=True, sharex=True, figsize=(6,6)):
    x = np.asarray(x)
    n_plots = len(y_sets)
    fig, axes = plt.subplots(n_plots, 1, sharex=sharex, figsize=figsize)
    if n_plots == 1:
        axes = [axes]

    # Colormap based on x
    cmap_obj = cm.get_cmap('viridis')
    norm = plt.Normalize(np.min(x), np.max(x))
    colors = cmap_obj(norm(x))

    for i, ax in enumerate(axes):
        y_group = y_sets[i]
        # Ensure we handle single arrays and lists of arrays
        if not isinstance(y_group[0], (list, np.ndarray)):
            y_group = [y_group]

        # Select label for this axis
        lbl = y_label[i] if isinstance(y_label, list) and i < len(y_label) else y_label

        # Plot each dataset in this axis
        for k, y in enumerate(y_group):
            y = np.array(y)
            label = None
            if isinstance(lbl, list):  # support nested label list
                label = lbl[k] if k < len(lbl) else None
            elif k == 0:
                label = lbl

            for j in range(len(x) - 1):
                if type_list is not None and type_list[j] == ['saddle', 'center', 'center']:
                    marker = 'o'
                    style  = '-'
                else:
                    marker = 's'
                    style  = '--'
                ax.plot(x[j:j+2], y[j:j+2], color=colors[j], linestyle=style, marker=marker)

            # Add dummy handle for legend
            if label:
                ax.plot([], [], color=cmap_obj(0.7), linestyle=style, label=label)

        # Grid and reference lines
        ax.grid(True, linestyle=":", linewidth=0.8)
        if hlines:
            for ypos, color, hlabel in hlines:
                ax.axhline(ypos, color=color, linestyle="--", label=hlabel)
        if vlines:
            for xpos, color, vlabel in vlines:
                ax.axvline(xpos, color=color, linestyle="--", label=vlabel)

        # Axis labels
        if isinstance(y_label, list) and not isinstance(y_label[i], list):
            ax.set_ylabel(y_label[i])
        elif isinstance(y_label, str):
            ax.set_ylabel(y_label)

        if legend:
            ax.legend()

    # Shared x-axis label and title
    if x_label:
        axes[-1].set_xlabel(x_label)
    if title:
        fig.suptitle(title)
    fig.tight_layout()
    return fig, axes

# Function to plot eigenvalues in the complex plane
def plot_eigenvalues_complex_plane(eigval_list, x0_list=None, pair_labels=None):
    n_pairs = len(eigval_list[0])//2
    fig, axes = plt.subplots(n_pairs, 1, constrained_layout=True)
    # Color map to encode x0 variation (or index)
    cmap = cm.viridis
    if x0_list is None:
        x0_list = np.arange(len(eigval_list))
    norm = plt.Normalize(np.min(x0_list), np.max(x0_list))
    colors = cmap(norm(x0_list))
    for j in range(n_pairs):
        ax = axes[j]
        ax.axhline(0, color='k', linewidth=0.8)
        ax.axvline(0, color='k', linewidth=0.8)
        ax.grid(True, linestyle=':', linewidth=0.8)
        for i, eigvals in enumerate(eigval_list):
            c = colors[i]
            ax.plot(eigvals[2*j].real, eigvals[2*j].imag, 'x', color=c, markersize=6)
            ax.plot(eigvals[2*j+1].real, eigvals[2*j+1].imag, 'x', color=c, markersize=6)
        ax.set_xlabel("Re($\lambda$)")
        ax.set_ylabel("Im($\lambda$)")
        label = pair_labels[j] if pair_labels else f"Pair {j+1}"
        ax.set_title(label)
        # Plot stability circle
        unit_circle = plt.Circle((0, 0), 1, facecolor='lightgrey', edgecolor='grey', linewidth=1)
        ax.add_patch(unit_circle)
    # Optional colorbar to show parameter variation (x0)
    if x0_list is not None:
        norm = plt.Normalize(min(x0_list), max(x0_list))
        sm = plt.cm.ScalarMappable(cmap='viridis', norm=norm)
        sm.set_array([])
        cbar = fig.colorbar(sm, ax=axes, orientation='horizontal', fraction=0.05, pad=0.1)
        cbar.set_label("$x_0$")
    return fig, axes

def jacobi_3d(state, mu):
    x = state[0,:]
    y = state[1,:]
    z = state[2,:]
    vx = state[3,:]
    vy = state[4,:]
    vz = state[5,:]
    r1 = ((x - mu)**2 + y**2 + z**2)**0.5
    r2 = ((x - mu + 1)**2 + y**2 + z**2)**0.5
    w = 0.5*(x**2+y**2) + (1-mu)/r1 + mu/r2
    return 2*w - (vx**2 + vy**2 + vz**2)

############
# RUN CODE #
############

x_ini = equilibrium_points['L1'][0]-0.002
x_end = equilibrium_points['P2'][0]+5e-5

# init_state_list = compute_solutions(x_ini, x_end, n=20, tol=1e-15, )
# np.savetxt("lyapunov.txt", init_state_list, delimiter=",")

init_state_list = load_solutions("lyapunov.txt", n=16)
full_state_list = integrate_solutions(init_state_list)

# Process results
x0_list     = []
vy0_list    = []
T_half_list = []
jacobi_list = []
for state in init_state_list:
    x0     = state[0]
    vy0    = state[4]
    T_half = state[6]
    x0_list.append(x0)
    vy0_list.append(vy0)
    T_half_list.append(T_half)
    jacobi_list.append(utils.jacobi_constant(x0, mu))

eigval_list = []
s1_list     = []
s2_list     = []
s3_list     = []
type_list   = []
for full_state in full_state_list:
    monodromy_matrix = full_state[6:,-1].reshape((6,6))
    s1, s2, s3, eigvals = stability_indices(monodromy_matrix)
    eigval_list.append(eigvals)
    # Use the stability index to set orbit type
    s_type = []
    for s in [s1,s2,s3]:
        if abs(s) <= 1.1:
            s_type.append('center')
        else:
            s_type.append('saddle')
    print(abs(s1),abs(s2),abs(s3),s_type)
    s1_list.append(abs(s1))
    s2_list.append(abs(s2))
    s3_list.append(abs(s3))
    type_list.append(s_type)

#########
# PLOTS #
#########

# i) Orbit plot
fig, ax = plt.subplots(1, 2, figsize=(8,4))
cmap = cm.viridis
x0_vals  = [s[0][0] for s in full_state_list]
norm = plt.Normalize(min(x0_vals), max(x0_vals))
for state, x0, s_type in zip(full_state_list, x0_vals, type_list):
    if s_type == ['saddle', 'center', 'center']:
        ls = '-'
    else:
        ls = '--'
    color = cmap(norm(x0))
    ax[0].plot(state[0], state[1], color=color, linestyle=ls, label=f"x0: {state[0][0]:.3g}, vy0: {state[4][0]:.3g}")
    ax[1].plot(state[3], state[4], color=color, linestyle=ls, label=f"x0: {state[0][0]:.3g}, vy0: {state[4][0]:.3g}")
# Equilibrium points
ax[0].plot(equilibrium_points['L1'][0], equilibrium_points['L1'][1], 'bX', label='L1')
ax[0].plot(equilibrium_points['L2'][0], equilibrium_points['L2'][1], 'rX', label='L2')
ax[0].plot(equilibrium_points['P2'][0], equilibrium_points['P2'][1], 'C1o', label='Jupiter')
# Axes formatting
for a in ax:
    a.axis("equal")
    a.legend()
    a.grid(True)
ax[0].set_xlabel("$x$")
ax[0].set_ylabel("$y$")
ax[1].set_xlabel("$\dot{x}$")
ax[1].set_ylabel("$\dot{y}$")

# Vertical lines for reference points
vlines = [(equilibrium_points['L1'][0], 'blue', 'L1'), (equilibrium_points['P2'][0], 'C1', 'Jupiter')]

# ii) x0 vs vy0 and x0 vs T/2
plot(x0_list, [vy0_list, T_half_list], y_label=["$\dot{y}_0$", "$T/2$"], x_label="$x_0$", vlines=vlines, type_list=type_list)

# x0 vs Re(λ) and Im(λ) for each pair
for j in range(3):
    eig_re_1 = [eigvals[2*j].real for eigvals in eigval_list]
    eig_im_1 = [eigvals[2*j].imag for eigvals in eigval_list]
    eig_re_2 = [eigvals[2*j+1].real for eigvals in eigval_list]
    eig_im_2 = [eigvals[2*j+1].imag for eigvals in eigval_list]
    plot(x0_list, [[eig_re_1, eig_re_2], [eig_im_1, eig_im_2]], y_label=["Re($\lambda$)", "Im($\lambda$)"], x_label="$x_0$", vlines=vlines, type_list=type_list, title=f"Eigenvalue pair {j+1}")

# x0 vs |s1|, |s2|, |s3|
plot(x0_list, [s1_list, s2_list, s3_list], y_label=["$|s_1|$", "$|s_2|$", "$|s_3|$"], x_label="$x_0$", hlines=[(1, "grey", "$|s_i|=1$")], vlines=vlines, type_list=type_list)

# x0 vs Cj and T
plot(x0_list, [jacobi_list, np.array(T_half_list)*2], y_label=["$C_j$", "$T$"], x_label="$x_0$", vlines=vlines, type_list=type_list)

# B) Eigenvalues on complex plane
plot_eigenvalues_complex_plane(eigval_list, x0_list=x0_list, pair_labels=[r"$(\lambda_1,\lambda_2)$", r"$(\lambda_3,\lambda_4)$", r"$(\lambda_5,\lambda_6)$"])

plt.show()
