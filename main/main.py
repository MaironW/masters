# Level 0 Module main
# Executes in order the Level 1 Modules DYN, SEN and NAV, storing information along the time
# Executes the PPC module to show results

from PPC           import PPC
from PPC.PPC_plots import PPC_plots
from SIM_par       import SIM_par
from Utils         import simulation

# Run simulation
timeline, DYN_obj, SEN_obj, NAV_obj = simulation.run()

# Stop simulation
simulation.stop()

# Plots
for key in SIM_par["PPC_plot_list"]:
    if key in PPC_plots:
        PPC_plots[key](timeline, DYN_obj, SEN_obj, NAV_obj)
PPC.show_plot()
