# Level 0 Module main
# Executes in order the Level 1 Modules DYN, SEN and NAV, storing information along the time
# Executes the PPC module to show results

from DYN   import DYN
from SEN   import SEN
from NAV   import NAV
from PPC   import PPC
from Utils import spice
from Utils import events
from Utils import integrator

# Simulation parameters
sim_dt         = 1      # [s]
sim_time_start = 0      # [s]
sim_time_end   = 3600*3 # [s] 3 h

sim_time       = sim_time_start # [s]
step           = 0
n_steps        = int((sim_time_end - sim_time_start)/sim_dt) + 1

# Initialize inputs table
events_table = events.build_events_table(sim_time_start, sim_time_end, sim_dt)

# Initialize outputs
DYN_out = DYN.initialize()
modules = DYN.get_dynamic_modules()

# Initialize timeline
timeline = PPC.init_timeline({"DYN" : DYN_out}, n_steps)

# Main loop
for step in range(n_steps):
    # Load external inputs
    inputs = events_table[sim_time]
    # Update algebraic modules first
    DYN_out = DYN.update_algebraic(sim_time, DYN_out, inputs["DYN"])
    # Integrate all dynamic states together
    DYN_out = integrator.rk4_step(sim_time, sim_dt, DYN_out, modules)
    # Save results into the timeline
    output = {"DYN" : DYN_out}
    PPC.update_timeline(timeline, output, step)

    # Update time
    sim_time += sim_dt

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
fig, ax = PPC.plot(timeline["DYN"]["DYN_TRA"]["SCpos_SSB"][:,0],      timeline["DYN"]["DYN_TRA"]["SCpos_SSB"][:,1],      timeline["DYN"]["DYN_TRA"]["SCpos_SSB"][:,2],      label="SC",    fig=fig, ax=ax)

# Plot inputs
fig, ax = PPC.plot(timeline["DYN"]["DYN_TIME"]["time_SIM"], timeline["DYN"]["DYN_GRV"]["grvacc_SSB"], label=["x","y","z"], xlabel="time_SIM [s]", ylabel="grvacc_SSB", title="grvacc_SSB")
fig, ax = PPC.plot(timeline["DYN"]["DYN_TIME"]["time_SIM"], timeline["DYN"]["DYN_TRA"]["SCpos_ECI"], label=["x","y","z"], xlabel="time_SIM [s]", ylabel="SCpos_ECI", title="SCpos_ECI")
fig, ax = PPC.plot(timeline["DYN"]["DYN_TIME"]["time_SIM"], timeline["DYN"]["DYN_TRA"]["SCpos_MCI"], label=["x","y","z"], xlabel="time_SIM [s]", ylabel="SCpos_MCI", title="SCpos_MCI")
fig, ax = PPC.plot(timeline["DYN"]["DYN_TIME"]["time_SIM"], timeline["DYN"]["DYN_TRA"]["SCpos_SCI"], label=["x","y","z"], xlabel="time_SIM [s]", ylabel="SCpos_SCI", title="SCpos_SCI")

fig, ax = PPC.plot(timeline["DYN"]["DYN_TRA"]["SCpos_TER"][:,0], timeline["DYN"]["DYN_TRA"]["SCpos_TER"][:,1], label=["x","y","z"], xlabel="SCpos_TER x [km]", ylabel="SCpos_TER y [km]", title="SCpos_TER")
fig, ax = PPC.plot(timeline["DYN"]["DYN_TRA"]["SCpos_MAR"][:,0], timeline["DYN"]["DYN_TRA"]["SCpos_MAR"][:,1], label=["x","y","z"], xlabel="SCpos_MAR x [km]", ylabel="SCpos_MAR y [km]", title="SCpos_MAR")

PPC.show_plot()
