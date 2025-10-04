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
def plot(x, y, z=None, style='', xlabel=None, ylabel=None, zlabel=None, label=None, title=None, fig=None, ax=None, subplot=None):
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

    # Plot data
    if z is None: # Plot 2D
        ax.plot(x, y, style, label=label)
        if xlabel: ax.set_xlabel(xlabel)
        if ylabel: ax.set_ylabel(ylabel)
    else: # Plot 3D
        ax.plot(x, y, z, style, label=label)
        if xlabel: ax.set_xlabel(xlabel)
        if ylabel: ax.set_ylabel(ylabel)
        if zlabel: ax.set_zlabel(zlabel)

    # Add legends
    handles, labels = ax.get_legend_handles_labels()
    if labels:
        ax.legend()

    if title:
        ax.set_title(title)
    ax.grid(True)

    return fig, ax

# Show plots after they are generated
def show_plot():
    plt.show()