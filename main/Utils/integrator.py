import numpy as np

def get_state_vector(state, modules):
    vec = []
    slices = {}
    i = 0
    for name, mod in modules.items():
        s = mod.get_state(state[name])
        n = len(s)
        slices[name] = slice(i, i+n)
        vec.append(s)
        i += n
    return np.hstack(vec), slices

def set_state_vector(state, y, slices, modules):
    for name, mod in modules.items():
        state[name] = mod.set_state(state[name], y[slices[name]])
    return state

# Perform one Runge-Kutta 4 step over all modules
def rk4_step(t, dt, state, modules):
    state_0, slices = get_state_vector(state, modules)

    # Set states, update algebraic modules and dynamic outputs, then evaluate derivatives
    def f(t_local, y_vec):
        state_tmp = set_state_vector(dict(state), y_vec, slices, modules)
        # Compute outputs of dynamic models
        for mod in modules.values():
            mod.outputs(t_local, state_tmp)
        # Compute derivatives of dynamic models
        dy = [mod.derivatives(t_local, state_tmp) for mod in modules.values()]
        return np.hstack(dy)

    k1 = f(t,        state_0)
    k2 = f(t + dt/2, state_0 + dt*k1/2)
    k3 = f(t + dt/2, state_0 + dt*k2/2)
    k4 = f(t + dt,   state_0 + dt*k3)
    state_next = state_0 + dt/6 * (k1 + 2*k2 + 2*k3 + k4)

    # Update global dictionary
    state = set_state_vector(state, state_next, slices, modules)

    return state
