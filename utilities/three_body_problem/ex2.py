import utils
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from scipy.integrate import solve_ivp

def variational_equations(t, state_M, mu):
    # Extract the state vector and the matrix M
    state  = state_M[:6]
    M_flat = state_M[6:]
    M      = M_flat.view().reshape(6,6)
    # Dynamic equations
    dstate = utils.state_equations(t, state, mu)
    # Jacobian matrix A(t)
    A = utils.linearization_matrix(state[0], state[1], state[2], mu)
    # Variational equations: dM/dt = A*M
    dM = A @ M
    dM_flat = dM.ravel()
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
    solution = solve_ivp(variational_equations, t_span, y0, args=(mu,), rtol=1e-12, atol=1e-12, events=y_event, dense_output=True)
    return solution

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
            state_half = solution.sol(t_end)[:6]
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
            M = solution.sol(t_end)[6:].reshape((6,6))
            dy_dvy  = M[1,4]
            dvx_dvy = M[3,4]

            # Build the minimal state transition matrix
            M_sub = np.array([[ dy_dvy, vy_half],
                              [dvx_dvy, ax_half]])

            # Correction of [vy0 only]
            delta = np.linalg.solve(M_sub, -F)
            state[4] += delta[0]

    return solution


# Mass of celestial bodies [kg]
body_mass = {
    "Sun"     : 1.988500e30,
    "Jupiter" : 1.898130e27,
}

# Compute the mass parameter µ for the Sun-Jupiter system
mu = utils.mu_from_masses(body_mass["Sun"], body_mass["Jupiter"])

# Compute equilibrium points
equilibrium_points = utils.equilibrium_points(mu)

# Get solutions
def compute_solutions():
    # Lists to store solutions
    init_state_list = []
    tf = 1.55
    vy = 0.01
    for x0 in np.linspace(equilibrium_points['L1'][0]-0.002, equilibrium_points['P2'][0]+5e-5, 20):
        # Brute force first guesses for vy0
        vy0_new = None
        for vy0 in np.linspace(vy, 0.6, 25):
            # vy0_new will be None until a valid solution is found.
            # In this case, keep using the brute force guesses.
            if vy0_new is not None:
                vy0 = vy0_new
            init_state = np.array([x0, 0, 0, 0, vy0, 0])
            solution = differential_correction(mu, init_state, tf)

            if solution.t_events[0].size > 0:
                T_half = solution.t_events[0][0]
                x_ini  = solution.sol(0)[0]
                x_half = solution.sol(T_half)[0]
                # If xf > x0, save updated init state
                if x_half > x_ini:
                    x,y,z,vx,vy,vz = solution.sol(0)[:6]
                    print(x,vy,T_half)
                    vy0_new = vy
                    init_state_list.append([x,y,z,vx,vy,vz,T_half])
                    tf = T_half*2
                    break # Stop iterating the velocity if solution found
        print(vy0_new)
    return init_state_list

# Load solutions from file in format [x,y,z,vx,vy,vz,T_half]
def load_solutions(file):
    init_state_list = np.loadtxt(file, delimiter=",")
    return init_state_list

# Integrate a full orbit given initial states
def integrate_solutions(init_state_list):
    full_state_list = []
    for init_state in init_state_list:
        x,y,z,vx,vy,vz,T_half = init_state
        tf = T_half
        solution = integrate_variational(mu, init_state[:6], 2*tf, terminate=False)
        t_vals = np.linspace(0, 2*tf, 2000)
        full_state = solution.sol(t_vals)
        full_state_list.append(full_state)
    return full_state_list

def stability_indices(M, tol=1e-6):
    M = np.asarray(M)
    eigvals = np.linalg.eigvals(M)
    eigvals = np.array(eigvals)
    # List pairs of eigvals
    used = np.zeros(len(eigvals), dtype=bool)
    pairs = []
    # Find reciprocs λ and 1/λ
    for i, lam in enumerate(eigvals):
        if used[i]:
            continue
        # Find λ' such that |λ * λ' - 1| < tol
        for j, lam2 in enumerate(eigvals):
            if i == j or used[j]:
                continue
            if np.abs(lam * lam2 - 1) < tol:
                pairs.append((lam, lam2))
                used[i] = True
                used[j] = True
                break
    # Sort the pairs per module to distinguish planar and transversal
    pair_mod = [np.mean([abs(lam[0]), abs(lam[1])]) for lam in pairs]
    idx_sorted = np.argsort(pair_mod)
    pairs_sorted = [pairs[i] for i in idx_sorted]
    # Assuming planar orbit
    planar_pairs = pairs_sorted[:2]
    transversal_pair = pairs_sorted[2]
    # Compute stability indices
    s1 = 0.5*(planar_pairs[0][0]  + 1/planar_pairs[0][0])
    s2 = 0.5*(planar_pairs[1][0]  + 1/planar_pairs[1][0])
    s3 = 0.5*(transversal_pair[0] + 1/transversal_pair[0])
    return s1, s2, s3, eigvals

# Generic plot function to show results
def plot(x, y_sets, x_label=None, y_label=None, vlines=None, hlines=None, title=None, markers=None, styles=None, legend=True, sharex=True, figsize=(6,6)):
    x = np.asarray(x)
    n_plots = len(y_sets)
    fig, axes = plt.subplots(n_plots, 1, sharex=sharex, figsize=figsize)
    if n_plots == 1:
        axes = [axes]
    # Create color mapping by x
    cmap_obj = cm.get_cmap('viridis')
    norm = plt.Normalize(np.min(x), np.max(x))
    colors = cmap_obj(norm(x))
    for i, ax in enumerate(axes):
        y = np.array(y_sets[i])
        lbl = y_label[i] if isinstance(y_label, list) and i < len(y_label) else y_label
        marker = markers[i] if markers and i < len(markers) else 'o'
        style = styles[i] if styles and i < len(styles) else '-'
        # Line with gradient color using small segments
        for j in range(len(x) - 1):
            ax.plot(x[j:j+2], y[j:j+2], color=colors[j], linestyle=style, marker=marker)
        if lbl:
            ax.plot([], [], color=cmap_obj(0.7), linestyle=style, label=lbl)
        # Grid and reference lines
        ax.grid(True, linestyle=":", linewidth=0.8)
        if hlines:
            for ypos, color, hlabel in hlines:
                ax.axhline(ypos, color=color, linestyle="--", label=hlabel)
        if vlines:
            for xpos, color, vlabel in vlines:
                ax.axvline(xpos, color=color, linestyle="--", label=vlabel)
        # Labels
        if isinstance(y_label, list):
            ax.set_ylabel(y_label[i])
        elif y_label:
            ax.set_ylabel(y_label)
        if legend and lbl:
            ax.legend()
    # Shared x-axis label
    if x_label:
        axes[-1].set_xlabel(x_label)
    if title:
        fig.suptitle(title)
    fig.tight_layout()
    return fig, axes

def plot_eigenvalues_complex_plane(eigval_list, x0_list=None, pair_labels=None):
    n_pairs = len(eigval_list[0]) // 2
    fig, axes = plt.subplots(1, n_pairs, figsize=(14, 4), constrained_layout=True)
    if n_pairs == 1:
        axes = [axes]
    # Color map to encode x0 variation (or index)
    cmap = cm.viridis
    colors = cmap(np.linspace(0, 1, len(eigval_list)))
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
    # Optional colorbar to show parameter variation (x0)
    if x0_list is not None:
        norm = plt.Normalize(min(x0_list), max(x0_list))
        sm = plt.cm.ScalarMappable(cmap='viridis', norm=norm)
        sm.set_array([])
        cbar = fig.colorbar(sm, ax=axes, orientation='horizontal', fraction=0.05, pad=0.1)
        cbar.set_label("$x_0$")
    return fig, axes

# init_state_list = compute_solutions()
# np.savetxt("lyapunov.txt", init_state_list, delimiter=",")

init_state_list = load_solutions("lyapunov.txt")
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
for full_state in full_state_list:
    monodromy_matrix = full_state.T[-1][6:].reshape((6,6))
    s1, s2, s3, eigvals = stability_indices(monodromy_matrix)
    eigval_list.append(eigvals)
    s1_list.append(abs(s1))
    s2_list.append(abs(s2))
    s3_list.append(abs(s3))

# i) Orbit plot
fig, ax = plt.subplots(1, 2, figsize=(8,4))
cmap = cm.viridis
x0_vals  = [s[0][0] for s in full_state_list]
norm = plt.Normalize(min(x0_vals), max(x0_vals))
for state, x0 in zip(full_state_list, x0_vals):
    color = cmap(norm(x0))
    ax[0].plot(state[0], state[1], color=color, label=f"x0: {state[0][0]:.3g}, vy0: {state[4][0]:.3g}")
    ax[1].plot(state[3], state[4], color=color, label=f"x0: {state[0][0]:.3g}, vy0: {state[4][0]:.3g}")
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
plot(x0_list, [vy0_list, T_half_list], y_label=["$\dot{y}_0$", "$T/2$"], x_label="$x_0$", vlines=vlines)

# x0 vs Re(λ) and Im(λ) for each pair
for j in range(3):
    eig_re = [eigvals[2*j].real for eigvals in eigval_list]
    eig_im = [eigvals[2*j].imag for eigvals in eigval_list]
    plot(x0_list, [eig_re, eig_im], y_label=["Re($\lambda$)", "Im($\lambda$)"], x_label="$x_0$", vlines=vlines, title=f"Eigenvalue pair {j+1}")

# x0 vs |s1|, |s2|, |s3|
plot(x0_list, [s1_list, s2_list, s3_list], y_label=["$|s_1|$", "$|s_2|$", "$|s_3|$"], x_label="$x_0$", hlines=[(1, "grey", "$|s_i|=1$")], vlines=vlines)

# x0 vs Cj and T
plot(x0_list, [jacobi_list, np.array(T_half_list)*2], y_label=["$C_j$", "$T$"], x_label="$x_0$", vlines=vlines)

# B) Eigenvalues on complex plane
plot_eigenvalues_complex_plane(eigval_list, x0_list=x0_list, pair_labels=[r"$(\lambda_1,\lambda_2)$", r"$(\lambda_3,\lambda_4)$", r"$(\lambda_5,\lambda_6)$"])

plt.show()
