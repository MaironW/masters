# Level 0 Module main
# Executes in order the Level 1 Modules DYN, SEN and NAV, storing information along the time
# Executes the PPC module to show results

from DYN.DYN import DYN
from PPC   import PPC
from Utils import spice
from Utils import events
from Utils import integrator

# Load data or run new simulation
DYN_log_save = False
DYN_log_load = False
DYN_log_path = "Logs/DYN"

# Simulation parameters
sim_dt         = 600            # [s] 10 min
sim_time_start = 0              # [s]
sim_time_end   = 3600*24*20     # [s] 20 days
sim_time       = sim_time_start # [s]
step           = 0
n_steps        = int((sim_time_end - sim_time_start)/sim_dt) + 1

# Initialize inputs table
events_table = events.build_events_table(sim_time_start, sim_time_end, sim_dt)

# Initialize Level-1 modules
DYN_obj = DYN()
DYN_out = DYN_obj.get_state
modules = DYN_obj.get_dynamic_modules()

# Initialize or load timeline
if DYN_log_load:
    DYN_timeline = PPC.load_timeline(DYN_log_path)
    timeline = {"DYN" : DYN_timeline}
else:
    timeline = {"DYN": PPC.init_timeline(DYN_obj.snapshot(), n_steps)}

# Main loop
for step in range(n_steps):
    # Load external inputs
    inputs = events_table[sim_time]

    # Load DYN from file or compute everything
    if DYN_log_load:
        DYN_out = PPC.load_module(timeline["DYN"], step)
    else:
        # Update algebraic modules first
        DYN_obj.update_algebraic(sim_time, inputs)
        # Integrate all dynamic states together
        DYN_out = integrator.rk4_step(sim_time, sim_dt, DYN_obj, inputs)
        # # Save results into the timeline
        output = {"DYN" : DYN_obj.snapshot()}
        PPC.update_timeline(timeline, output, step)

    # Update time
    sim_time += sim_dt

# Save log (currently only for DYN module)
if DYN_log_save:
    PPC.store_timeline(timeline, DYN_log_path)

# Clear Kernels from memory
spice.clear_kernels()

# # Plot trajectory
fig, ax = PPC.plot([0], [0], [0], style='+', label="SSB", xlabel="X SSB [km]", ylabel="Y SSB [km]", zlabel="Z SSB [km]", title="Simulation Trajectories")
fig, ax = PPC.plot(timeline["DYN"]["DYN_SUN"]["SUNpos_SSB"][:,0],     timeline["DYN"]["DYN_SUN"]["SUNpos_SSB"][:,1],     timeline["DYN"]["DYN_SUN"]["SUNpos_SSB"][:,2],     label="Sun",   fig=fig, ax=ax)
fig, ax = PPC.plot(timeline["DYN"]["DYN_EARTH"]["EARTHpos_SSB"][:,0], timeline["DYN"]["DYN_EARTH"]["EARTHpos_SSB"][:,1], timeline["DYN"]["DYN_EARTH"]["EARTHpos_SSB"][:,2], label="Earth", fig=fig, ax=ax)
fig, ax = PPC.plot(timeline["DYN"]["DYN_EARTH"]["MOONpos_SSB"][:,0],  timeline["DYN"]["DYN_EARTH"]["MOONpos_SSB"][:,1],  timeline["DYN"]["DYN_EARTH"]["MOONpos_SSB"][:,2],  label="Moon",  fig=fig, ax=ax)
fig, ax = PPC.plot(timeline["DYN"]["DYN_MARS"]["MARSpos_SSB"][:,0],   timeline["DYN"]["DYN_MARS"]["MARSpos_SSB"][:,1],   timeline["DYN"]["DYN_MARS"]["MARSpos_SSB"][:,2],   label="Mars",  fig=fig, ax=ax)
fig, ax = PPC.plot(timeline["DYN"]["DYN_MARS"]["DEIMOSpos_SSB"][:,0], timeline["DYN"]["DYN_MARS"]["DEIMOSpos_SSB"][:,1], timeline["DYN"]["DYN_MARS"]["DEIMOSpos_SSB"][:,2], label="Deimos",fig=fig, ax=ax)
fig, ax = PPC.plot(timeline["DYN"]["DYN_MARS"]["PHOBOSpos_SSB"][:,0], timeline["DYN"]["DYN_MARS"]["PHOBOSpos_SSB"][:,1], timeline["DYN"]["DYN_MARS"]["PHOBOSpos_SSB"][:,2], label="Phobos",fig=fig, ax=ax)
fig, ax = PPC.plot(timeline["DYN"]["DYN_TRA"]["SCpos_SSB"][:,0],      timeline["DYN"]["DYN_TRA"]["SCpos_SSB"][:,1],      timeline["DYN"]["DYN_TRA"]["SCpos_SSB"][:,2],      label="SC",    fig=fig, ax=ax)

# # Plot 2D
fig, ax = PPC.plot(timeline["DYN"]["DYN_TIME"]["time_SIM"], timeline["DYN"]["DYN_TIME"]["time_UTC"], xlabel="time_SIM [s]", ylabel="time_UTC [s]", title="time")
fig, ax = PPC.plot(timeline["DYN"]["DYN_TIME"]["time_SIM"], timeline["DYN"]["DYN_SUN"]["SUNpos_SSB"], xlabel="time_SIM [s]", ylabel="SUNpos_SSB [km]", title="SUNpos")
fig, ax = PPC.plot(timeline["DYN"]["DYN_TIME"]["time_SIM"], timeline["DYN"]["DYN_EARTH"]["EARTHpos_SSB"], xlabel="time_SIM [s]", ylabel="EARTHpos_SSB [km]", title="EARTHpos")
fig, ax = PPC.plot(timeline["DYN"]["DYN_TIME"]["time_SIM"], timeline["DYN"]["DYN_EARTH"]["MOONpos_SSB"], xlabel="time_SIM [s]", ylabel="MOONpos_SSB [km]", title="MOONpos", fig=fig, ax=ax)
fig, ax = PPC.plot(timeline["DYN"]["DYN_TIME"]["time_SIM"], timeline["DYN"]["DYN_MARS"]["MARSpos_SSB"], xlabel="time_SIM [s]", ylabel="MARSpos_SSB [km]", title="MARSpos")
fig, ax = PPC.plot(timeline["DYN"]["DYN_TIME"]["time_SIM"], timeline["DYN"]["DYN_ATT"]["SSBq_BOF"], xlabel="time_SIM [s]", ylabel="SSBq_BOF", title="SSBq_BOF")
fig, ax = PPC.plot(timeline["DYN"]["DYN_TIME"]["time_SIM"], timeline["DYN"]["DYN_TRA"]["SCpos_ECI"], xlabel="time_SIM [s]", ylabel="SCpos_ECI [km]", title="SCpos_ECI")
fig, ax = PPC.plot(timeline["DYN"]["DYN_TIME"]["time_SIM"], timeline["DYN"]["DYN_GRV"]["grvacc_SSB"], xlabel="time_SIM [s]", ylabel="grvacc_SSB [km/s^2]", title="grvacc_SSB")

PPC.show_plot()
