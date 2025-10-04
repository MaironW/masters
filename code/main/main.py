# Level 0 Module main
# Executes in order the Level 1 Modules DYN, SEN and NAV, storing information along the time
# Executes the PPC module to show results

from DYN   import DYN
from SEN   import SEN
from NAV   import NAV
from PPC   import PPC
from Utils import spice

# Simulation parameters
sim_dt         = 3600*24     # [s] 1 day
sim_time_start = 0           # [s]
sim_time_end   = 3600*24*365 # [s] 365 days
step           = 0
n_steps        = int((sim_time_end - sim_time_start)/sim_dt) + 1

# Initialize DYN variables
DYN_out = DYN.DYN_out

# Initialize timeline
timeline = PPC.init_timeline({"DYN" : DYN.DYN_out, "SEN" : SEN.SEN_out, "NAV" : NAV.NAV_out}, n_steps)

# Main loop
for step in range(n_steps):
    DYN_out = DYN.run()
    SEN_out = SEN.run(DYN_out)
    NAV_out = NAV.run(SEN_out)

    # Save results into the timeline
    output = {"DYN" : DYN_out, "SEN" : SEN_out, "NAV" : NAV_out}
    PPC.update_timeline(timeline, output, step)

# Clear Kernels from memory
spice.clear_kernels()

# Plot trajectory
fig, ax = PPC.plot([0], [0], [0], style='+', label="SSB", xlabel="X SSB [km]", ylabel="Y SSB [km]", zlabel="Z SSB [km]", title="Simulation Trajectories")
fig, ax = PPC.plot(timeline["DYN"]["DYN_SUN"]["SUNpos_SSB"][:,0],     timeline["DYN"]["DYN_SUN"]["SUNpos_SSB"][:,1],     timeline["DYN"]["DYN_SUN"]["SUNpos_SSB"][:,2],     label="Sun",   fig=fig, ax=ax)
fig, ax = PPC.plot(timeline["DYN"]["DYN_EARTH"]["EARTHpos_SSB"][:,0], timeline["DYN"]["DYN_EARTH"]["EARTHpos_SSB"][:,1], timeline["DYN"]["DYN_EARTH"]["EARTHpos_SSB"][:,2], label="Earth", fig=fig, ax=ax)
fig, ax = PPC.plot(timeline["DYN"]["DYN_EARTH"]["MOONpos_SSB"][:,0],  timeline["DYN"]["DYN_EARTH"]["MOONpos_SSB"][:,1],  timeline["DYN"]["DYN_EARTH"]["MOONpos_SSB"][:,2],  label="Moon",  fig=fig, ax=ax)
fig, ax = PPC.plot(timeline["DYN"]["DYN_MARS"]["MARSpos_SSB"][:,0],   timeline["DYN"]["DYN_MARS"]["MARSpos_SSB"][:,1],   timeline["DYN"]["DYN_MARS"]["MARSpos_SSB"][:,2],   label="Mars",  fig=fig, ax=ax)
fig, ax = PPC.plot(timeline["DYN"]["DYN_MARS"]["DEIMOSpos_SSB"][:,0], timeline["DYN"]["DYN_MARS"]["DEIMOSpos_SSB"][:,1], timeline["DYN"]["DYN_MARS"]["DEIMOSpos_SSB"][:,2], label="Deimos",fig=fig, ax=ax)
fig, ax = PPC.plot(timeline["DYN"]["DYN_MARS"]["PHOBOSpos_SSB"][:,0], timeline["DYN"]["DYN_MARS"]["PHOBOSpos_SSB"][:,1], timeline["DYN"]["DYN_MARS"]["PHOBOSpos_SSB"][:,2], label="Phobos",fig=fig, ax=ax)

PPC.show_plot()
