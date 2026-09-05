# Module to organize outputs and show results

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import os

# Define color cycle
colors = {
    "blue"      : "#0066FF",
    "red"       : "#CC0000",
    "green"     : "#33CC00",
    "magenta"   : "#FF00FF",
    "orange"    : "#FE9920",
    "purple"    : "#8052CF",
    "yellow"    : "#F5F22B",
    "grey"      : "#505050",
    "lightgrey" : "#707070",
    "darkgrey"  : "#303030",
    "black"     : "#111111",
}

# Initialize timeline to store data for post processing
# Scalars  -> (n_steps,)
# Vectors  -> (n_steps, dim)
# Matrices -> (n_steps, rows, cols)
def init_timeline(sample, n_steps):
    def recurse(node):
        if isinstance(node, dict):
            return {k: recurse(v) for k, v in node.items()}
        arr = np.array(node)
        # Case numeric
        if np.issubdtype(arr.dtype, np.number):
            out = np.full((n_steps, *arr.shape), np.nan)
            out[0] = arr
            return out
        # Case non-numeric (strings, objects, etc.)
        else:
            out = np.empty((n_steps, *arr.shape), dtype=object)
            out[:] = None
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
    os.makedirs(path, exist_ok=True)

    # Get highest module to store
    modules_order = ["DYN", "SEN", "NAV"]

    # Detect level automatically from timeline keys
    present = [m for m in modules_order if m in timeline]
    if not present:
        return

    max_idx = max(modules_order.index(m) for m in present)
    data_to_save = {}

    # Save up to highest level
    for m in modules_order[:max_idx + 1]:
        if m not in timeline:
            continue
        data_to_save[m] = np.array(timeline[m], dtype=object)
    np.savez_compressed(path+"/timeline.npz", **data_to_save)

# Load timeline data from file
def load_timeline(timeline, modules, path):
    if not os.path.exists(path):
        return timeline
    if not modules:
        return timeline

    data = np.load(path, allow_pickle=True)

    # Load up to highest level
    for m in modules:
        if m not in data:
            continue
        timeline[m] = data[m].item()
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
    # Create a new figure
    if fig is None:
        fig = plt.figure()

    # Create or select subplot/axes
    if ax is None:
        if subplot is not None:
            nrows, ncols, index = subplot

            # Initialize shared axis reference if not present
            if not hasattr(fig, "_shared_x_axis"):
                fig._shared_x_axis = None
            if z is None:
                if fig._shared_x_axis is None:
                    ax = fig.add_subplot(nrows, ncols, index)
                    fig._shared_x_axis = ax  # first subplot becomes reference
                else:
                    ax = fig.add_subplot(nrows, ncols, index, sharex=fig._shared_x_axis)
            else:
                if fig._shared_x_axis is None:
                    ax = fig.add_subplot(nrows, ncols, index, projection='3d')
                    fig._shared_x_axis = ax
                else:
                    ax = fig.add_subplot(nrows, ncols, index, projection='3d', sharex=fig._shared_x_axis)
        else:
            if z is None:
                ax = fig.gca()
            else:
                ax = fig.add_subplot(111, projection='3d')

    # Apply colors
    kwargs = {}
    if color is not None:
        kwargs['color'] = color

    # Apply zorder
    if zorder is not None:
        kwargs['zorder'] = zorder

    # Convert scalars to 1-element arrays
    if np.isscalar(x):
        x = np.atleast_1d(x)

    if np.isscalar(y):
        y = np.atleast_1d(y)

    if z is not None and np.isscalar(z):
        z = np.atleast_1d(z)

    # Plot data
    if z is None: # Plot 2D
        if len(x) > 0: # Only plot when there is data, to keep the color sequence
            ax.plot(x, y, style, **kwargs, label=label)
        if xlabel: ax.set_xlabel(xlabel)
        if ylabel: ax.set_ylabel(ylabel)
    else: # Plot 3D
        if len(x) > 0: # Only plot when there is data, to keep the color sequence
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

# Function to allow saving multiple objects into a json file
# Useful for storing parameters of Monte Carlo runs
def serialize_json(obj):
    if isinstance(obj, np.ndarray):
        return obj.tolist()
    if isinstance(obj, np.generic):  # numpy scalar
        return obj.item()
    if isinstance(obj, dict):
        return {k: serialize_json(v) for k, v in obj.items()}
    if isinstance(obj, list):
        return [serialize_json(v) for v in obj]
    return obj

# Downsample array of data from the timeline
# Useful to generate tikz figures
def downsample(original_array, n_points):
    idx = np.linspace(0, len(original_array)-1, n_points, dtype=int)
    array_downsampled = original_array[idx]
    return array_downsampled