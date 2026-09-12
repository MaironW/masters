import copy
import numpy as np

# Concatenate list of vectors into one
def _flatten(vecs):
    if not vecs:
        return np.array([])
    return np.concatenate(vecs)

# Split flat vector into list of vectors with given shapes
def _unflatten(vec, shapes):
    out = []
    i = 0
    for n in shapes:
        out.append(vec[i:i+n])
        i += n
    return out

# Perform one Runge-Kutta 4 step over all modules
def rk4_step(t, dt, state, inputs):
    # Pack dynamic state
    state_0_list = state.get_state()
    sizes = [len(v) for v in state_0_list]
    state_0 = _flatten(state_0_list)

    # Create a temporary copy of the state just to operate inside f()
    # without changing the final state
    # It is fine to do it only once because the function f() is returning a
    # flattened version of the new states, which is by itself a copy
    state_tmp = copy.copy(state)

    # Set states, update algebraic modules and dynamic outputs, then evaluate derivatives
    def f(t_local, y_vec):
        vecs = _unflatten(y_vec, sizes)

        # Inject dynamic state
        state_tmp.set_state(vecs)
        state_tmp.update_algebraic(t_local, None, inputs)

        # Compute derivatives of dynamic models
        dy = state_tmp.derivatives(t_local)
        return _flatten(dy)

    k1 = f(t,        state_0)
    k2 = f(t + dt/2, state_0 + dt*k1/2)
    k3 = f(t + dt/2, state_0 + dt*k2/2)
    k4 = f(t + dt,   state_0 + dt*k3)
    state_next = state_0 + dt/6 * (k1 + 2*k2 + 2*k3 + k4)

    # Write back final state
    state.set_state(_unflatten(state_next, sizes))

    return state
