# Module to organize outputs and show results

import numpy as np
import matplotlib.pyplot as plt

# Initialize timeline to store data for post processing
# Scalars  -> (n_steps,)
# Vectors  -> (n_steps, dim)
# Matrices -> (n_steps, rows, cols)
def init_timeline(sample, n_steps):
    def recurse(node):
        if isinstance(node, dict):
            return {k: recurse(v) for k, v in node.items()}
        else:
            arr = np.array(node)
            return np.full((n_steps, *arr.shape), np.nan)
    return recurse(sample)

# Fill preallocated timeline arrays with values at index=step.
def update_timeline(timeline, output, step):
    def recurse(t_node, o_node, idx):
        if isinstance(o_node, dict):
            for k, v in o_node.items():
                if k in t_node:
                    recurse(t_node[k], v, idx)
                # if missing in output, leave NaN in place
        else:
            t_node[idx] = o_node
    recurse(timeline, output, step)

# Plot function
def plot(x, y, xlabel=None, ylabel=None, title=None):
    plt.plot(x, y)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)
    plt.grid(True)
    plt.show()
