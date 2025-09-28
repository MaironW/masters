# Level 0 Module main
# Executes in order the Level 1 Modules DYN, SEN and NAV, storing information along the time
# Executes the PPC module to show results

from DYN import DYN
from SEN import SEN
from NAV import NAV
from PPC.PPC import plot_series

# Simulation parameters
sim_dt         = 1 # [s]
sim_time_start = 0   # [s]
sim_time_end   = 1000 # [s]
sim_time       = sim_time_start
iter           = 0
n_iters        = int((sim_time_end - sim_time_start)/sim_dt) + 1
timeline       = [0]*n_iters

# Initialize DYN variables
DYN_out = DYN.DYN_out

# Main loop
while sim_time < sim_time_end:
    DYN_out = DYN.run()
    SEN_out = SEN.run(DYN_out)
    NAV_out = NAV.run(SEN_out)

    timeline[iter] = {"DYN" : DYN_out, "SEN" : SEN_out, "NAV" : NAV_out}
    # Move on time
    iter     += 1
    sim_time += sim_dt # [s]

plot_series(timeline, "DYN.DYN_TIME.time", "SEN.SEN_STR.time")