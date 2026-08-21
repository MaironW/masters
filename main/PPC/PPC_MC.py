# Monte Carlo PPC Module

from PPC import PPC
import numpy as np
import glob
import os
import matplotlib.pyplot as plt

# Extract a variable from a nested timeline dictionary.
def _get_variable(data, path, axis=None):
    if isinstance(path, str):
        path = path.split(".")
    node = data
    for key in path:
        if not isinstance(node, dict):
            raise KeyError(f"Cannot access '{key}': current object is not a dictionary.")
        if key not in node:
            raise KeyError(f"Variable '{key}' not found. Available keys: {list(node.keys())}")
        node = node[key]
    arr = np.asarray(node)

    # If axis/component index is provided, slice along the last dimension
    if axis is not None:
        arr = arr[..., axis]

    return arr

# Recursively load all Monte Carlo .npz files below `path`.
def load_mc_logs(path, pattern="*.npz"):
    search_path = os.path.join(path, "**", pattern)
    files = sorted(glob.glob(search_path, recursive=True))

    if len(files) == 0:
        raise FileNotFoundError(f"No '{pattern}' files found below:\n{path}")

    runs = []
    for filename in files:
        data = np.load(filename, allow_pickle=True)
        timeline = {}
        for key in data.files:
            timeline[key] = data[key].item()
        runs.append(timeline)
        data.close()

    print(f"Loaded {len(runs)} Monte Carlo runs.")
    return runs, files


# Extract x/y data from all Monte Carlo runs and post process y data
def mc_data(runs, postprocess, x_var, x_axis=None, y_axis=None, check_x=True):
    x_all = []
    y_all = []

    for i, timeline in enumerate(runs):
        x = _get_variable(timeline, x_var, axis=x_axis)
        y = postprocess(timeline)
        x = np.asarray(x)
        y = np.asarray(y)

        if y_axis is not None:
            y = y[..., y_axis]

        # Remove singleton dimensions
        x = np.squeeze(x)
        y = np.squeeze(y)

        if x.ndim != 1:
            raise ValueError(f"x variable '{x_var}' in run {i} is not 1D. Shape = {x.shape}")

        # Allow vector-valued variables only if they are explicitly selected beforehand.
        if y.ndim != 1:
            raise ValueError(f"Processed y variable in run {i} is not 1D. Shape = {y.shape}\n")

        if len(x) != len(y):
            raise ValueError(f"x and y have different lengths in run {i}: {len(x)} vs {len(y)}")

        x_all.append(x)
        y_all.append(y)

    # Check that all runs use the same x axis
    if check_x:
        for i in range(1, len(x_all)):
            if not np.allclose(x_all[0], x_all[i], equal_nan=True):
                raise ValueError(f"x variable '{x_var}' is different between run 0 and run {i}.")

    x = x_all[0]

    # Stack runs
    y_runs = np.vstack(y_all)

    # Average independently at each timestep.
    # NaNs are ignored.
    y_mean = np.nanmean(y_runs, axis=0)
    y_std = np.nanstd(y_runs, axis=0)

    return x, y_runs, y_mean, y_std

# Plot all individual Monte Carlo runs together
def plot_mc(x_var, y_runs,
            x_axis=None, y_axis=None,
            xlabel=None, ylabel=None,
            title=None, color=None, alpha=0.15,
            run_label=None, fig=None, ax=None, subplot=None):

    # Defaults
    if color is None:
        color = PPC.colors["blue"]

    if xlabel is None:
        xlabel = f"{x_var}[{x_axis}]" if x_axis is not None else x_var

    if ylabel is None:
        ylabel = f"{y_runs}[{y_axis}]" if y_axis is not None else y_runs

    if title is None:
        title = f"Monte Carlo: {y_runs}"

    if run_label is None:
        run_label = f"Monte Carlo runs"

    # Create figure
    fig, ax = PPC.plot([], [],  xlabel=xlabel, ylabel=ylabel, title=title, fig=fig, ax=ax, subplot=subplot)

    # Plot individual Monte Carlo runs
    for i in range(y_runs.shape[0]):
        plot_label = run_label if i == 0 else None
        fig, ax = PPC.plot(x_var, y_runs[i], color=color, label=plot_label, fig=fig, ax=ax)

        # Modify the resulting Line2D object to make the Monte Carlo realizations transparent.
        line = ax.lines[-1]
        line.set_alpha(alpha)

    # Rebuild legend after modifying the lines
    ax.legend()

    return fig, ax
