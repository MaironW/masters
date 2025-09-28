# Module to organize outputs and show results

import matplotlib.pyplot as plt

# Extract a list of values from a nested dict structure given a dot-separated path.
# Example: "DYN.DYN_TIME.time"
def extract_series(timeline, path):
    keys = path.split(".")
    series = []
    for snap in timeline:
        val = snap
        try:
            for k in keys:
                val = val[k]
        except (KeyError, TypeError):
            val = None   # pad with None if missing
        series.append(val)
    return series

def plot_series(timeline, x_path, y_path, xlabel=None, ylabel=None, title=None):
    x = extract_series(timeline, x_path)
    y = extract_series(timeline, y_path)

    plt.plot(x, y)
    plt.xlabel(xlabel or x_path)
    plt.ylabel(ylabel or y_path)
    plt.title(title or f"{y_path} vs {x_path}")
    plt.grid(True)
    plt.show()