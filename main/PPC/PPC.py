# Module to organize outputs and show results

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import os

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
            out = np.full((n_steps, *arr.shape), np.nan)
            out[0] = arr
            return out
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

# Store the timeline data to be loaded in the future
def store_timeline(timeline, path):
    for module_name, module_dict in timeline["DYN"].items():
        os.makedirs(path, exist_ok=True)
        filename = os.path.join(path, f"{module_name}.npz")
        np.savez_compressed(filename, data=module_dict)

# Load timeline data from file
def load_timeline(path):
    timeline = {}
    for fname in os.listdir(path):
        if not fname.endswith(".npz"):
            continue
        module_name = fname[:-4] # remove ".npz"
        full_path = os.path.join(path, fname)
        data = np.load(full_path, allow_pickle=True)["data"].item()
        timeline[module_name] = data
    return timeline

# Extracts one module (DYN, SEN, NAV) for a given step,
# preserving the hierarchical structure exactly as stored.
def load_module(timeline, step):
    def recurse(node):
        # Case 1: nested dict → recurse
        if isinstance(node, dict):
            return {k: recurse(v) for k, v in node.items()}
        # Case 2: numpy array → extract the step-th entry
        elif isinstance(node, np.ndarray):
            return node[step]
        # Safety: anything else is returned as-is
        else:
            return node
    return recurse(timeline)

# Change default matplotlib configuration
def setup_plot(colors):
    mpl.rcParams['axes.prop_cycle'] = mpl.cycler(color=colors)

# Plot function
def plot(x, y, z=None, style='', color=None, xlabel=None, ylabel=None, zlabel=None, label=None, title=None, fig=None, ax=None, subplot=None, aspect=None, zorder=None):
    if fig is None:
        fig = plt.figure()

    # Create or select subplot/axes
    if ax is None:
        if subplot is not None:
            nrows, ncols, index = subplot
            if z is None:
                ax = fig.add_subplot(nrows, ncols, index)
            else:
                ax = fig.add_subplot(nrows, ncols, index, projection='3d')
        else:
            if z is None:
                ax = fig.gca() # get current axes (2D)
            else:
                ax = fig.add_subplot(111, projection='3d')

    # Apply colors
    kwargs = {}
    if color is not None:
        kwargs['color'] = color

    # Apply zorder
    if zorder is not None:
        kwargs['zorder'] = zorder

    # Plot data
    if z is None: # Plot 2D
        ax.plot(x, y, style, **kwargs, label=label)
        if xlabel: ax.set_xlabel(xlabel)
        if ylabel: ax.set_ylabel(ylabel)
    else: # Plot 3D
        ax.plot(x, y, z, style, **kwargs, label=label)
        if xlabel: ax.set_xlabel(xlabel)
        if ylabel: ax.set_ylabel(ylabel)
        if zlabel: ax.set_zlabel(zlabel)

    # Add legends
    handles, labels = ax.get_legend_handles_labels()
    if labels:
        ax.legend()

    if title:
        ax.set_title(title)

    if aspect:
        ax.set_aspect(aspect)

    ax.grid(True)

    return fig, ax

# Show plots after they are generated
def show_plot():
    plt.show()

# Generate a gnomonic projection of a vector
def gnomonic_projection(vec):
    x, y, z = vec
    mask = z > 0
    u = np.full_like(z, np.nan)
    v = np.full_like(z, np.nan)
    u[mask] = x[mask] / z[mask]
    v[mask] = y[mask] / z[mask]
    return u, v
