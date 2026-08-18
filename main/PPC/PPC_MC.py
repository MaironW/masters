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


# Extract x/y data from all Monte Carlo runs
# and calculate the mean y value at every simulation step.
def mc_data(runs, x_var, y_var, x_axis=None, y_axis=None, check_x=True):
    x_all = []
    y_all = []

    for i, run in enumerate(runs):
        x = _get_variable(run, x_var, axis=x_axis)
        y = _get_variable(run, y_var, axis=y_axis)
        x = np.asarray(x)
        y = np.asarray(y)

        # Remove singleton dimensions
        x = np.squeeze(x)
        y = np.squeeze(y)

        if x.ndim != 1:
            raise ValueError(f"x variable '{x_var}' in run {i} is not 1D. Shape = {x.shape}")

        # Allow vector-valued variables only if they are explicitly
        # selected beforehand.
        if y.ndim != 1:
            raise ValueError(f"y variable '{y_var}' in run {i} is not 1D. Shape = {y.shape}\n"
                "Select one component of the variable first.")

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

    return x, y_runs, y_mean

# Load Monte Carlo logs and plot all individual runs together with their mean.
def plot_mc(log_path, 
            x_var, y_var, 
            x_axis=None, y_axis=None, 
            xlabel=None, ylabel=None, 
            title=None, color_runs=None, color_mean=None, alpha=0.15, mean_linewidth=2.5, 
            run_label=None):
    
    # Defaults
    if color_runs is None:
        color_runs = PPC.colors["blue"]

    if color_mean is None:
        color_mean = PPC.colors["blue"]

    if xlabel is None:
        xlabel = f"{x_var}[{x_axis}]" if x_axis is not None else x_var

    if ylabel is None:
        ylabel = f"{y_var}[{y_axis}]" if y_axis is not None else y_var

    if title is None:
        title = f"Monte Carlo: {y_var}"

    if run_label is None:
        run_label = f"y(t) {os.path.basename(os.path.normpath(log_path))}"

    # Load data
    runs, files = load_mc_logs(log_path)

    x, y_runs, y_mean = mc_data(
        runs,
        x_var=x_var,
        y_var=y_var,
        x_axis=x_axis,
        y_axis=y_axis
    )

    # Create figure
    # THIS COULD BE DONE BY THE PPC.Plot
    fig = plt.figure()
    ax = None

    # Plot individual Monte Carlo runs
    for i in range(y_runs.shape[0]):
        plot_label = run_label if i == 0 else None
        fig, ax = PPC.plot(x, y_runs[i], color=color_runs, label=plot_label, fig=fig, ax=ax)

        # Modify the resulting Line2D object to make the Monte Carlo realizations transparent.
        line = ax.lines[-1]
        line.set_alpha(alpha)

    # Plot mean
    fig, ax = PPC.plot(x, y_mean, color=color_mean, label="y(t) Real Mean", fig=fig, ax=ax)

    mean_line = ax.lines[-1]
    mean_line.set_linewidth(mean_linewidth)

    # Labels
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title)

    # Rebuild legend after modifying the lines
    ax.legend()

    return fig, ax, x, y_runs, y_mean
