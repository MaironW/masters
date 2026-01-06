# Level 0 Module main
# Executes in order the Level 1 Modules DYN, SEN and NAV, storing information along the time
# Executes the PPC module to show results

from DYN.DYN       import DYN
from SEN.SEN       import SEN
from PPC           import PPC
from PPC.PPC_plots import PPC_plots
from SIM_par       import SIM_par
from Utils         import spice
from Utils         import events
from Utils         import integrator

# Load data or run new simulation
DYN_log_save = SIM_par["DYN_log_save"]
DYN_log_load = SIM_par["DYN_log_load"]
DYN_log_path = SIM_par["DYN_log_path"]

# Simulation parameters
sim_dt         = SIM_par["dt"]
sim_time_start = SIM_par["time_start"]
sim_time_end   = SIM_par["time_end"]
sim_time       = sim_time_start
step           = 0
n_steps        = int((sim_time_end - sim_time_start)/sim_dt) + 1

# Initialize inputs table
events_table = events.build_events_table(sim_time_start, sim_time_end, sim_dt)

# Initialize Level-1 modules
DYN_obj = DYN()
SEN_obj = SEN()

# Initialize or load timeline
if DYN_log_load:
    DYN_timeline = PPC.load_timeline(DYN_log_path)
    timeline = {
        "DYN" : DYN_timeline,
        "SEN": PPC.init_timeline(SEN_obj.snapshot(), n_steps),
    }
else:
    timeline = {
        "DYN": PPC.init_timeline(DYN_obj.snapshot(), n_steps),
        "SEN": PPC.init_timeline(SEN_obj.snapshot(), n_steps),
    }

# Main loop
for step in range(1, n_steps):
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
        # Update SEN
        SEN_obj.update_algebraic(sim_time, inputs)
        # Save results into the timeline
        output = {
            "DYN" : DYN_obj.snapshot(),
            "SEN" : SEN_obj.snapshot()
        }
        PPC.update_timeline(timeline, output, step)

    # Update time
    sim_time += sim_dt

# Save log (currently only for DYN module)
if DYN_log_save:
    PPC.store_timeline(timeline, DYN_log_path)

# Clear Kernels from memory
spice.clear_kernels()

# Plots
for key in SIM_par["PPC_plot_list"]:
    if key in PPC_plots:
        PPC_plots[key](timeline)
PPC.show_plot()
