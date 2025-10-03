# Level 0 Module main
# Executes in order the Level 1 Modules DYN, SEN and NAV, storing information along the time
# Executes the PPC module to show results

from DYN import DYN
from SEN import SEN
from NAV import NAV
from PPC import PPC

# Simulation parameters
sim_dt         = 1    # [s]
sim_time_start = 0    # [s]
sim_time_end   = 2    # [s]
sim_time       = sim_time_start
step           = 0
n_steps        = int((sim_time_end - sim_time_start)/sim_dt) + 1

# Initialize timeline
timeline = PPC.init_timeline({"DYN" : DYN.DYN_out, "SEN" : SEN.SEN_out, "NAV" : NAV.NAV_out}, n_steps)

# Initialize DYN variables
DYN_out = DYN.DYN_out

# Main loop
while sim_time <= sim_time_end:
    DYN_out = DYN.run()
    SEN_out = SEN.run(DYN_out)
    NAV_out = NAV.run(SEN_out)

    # Save results into the timeline
    output = {"DYN" : DYN_out, "SEN" : SEN_out, "NAV" : NAV_out}
    PPC.update_timeline(timeline, output, step)

    # Move on time
    step     += 1
    sim_time += sim_dt # [s]

print(timeline)