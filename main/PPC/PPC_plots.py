# List all plots to be generated, organized by group name, to be selected on SIM_par

from PPC import PPC

from .plots.DYN_plots import DYN_plots
from .plots.SEN_plots import SEN_plots
from .plots.NAV_plots import NAV_plots

PPC.setup_plot(PPC.colors.values())

PPC_plots = DYN_plots | SEN_plots | NAV_plots
