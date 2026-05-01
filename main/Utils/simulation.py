# Utils functions to run the simulations

import time
from DYN.DYN import DYN
from SEN.SEN import SEN
from NAV.NAV import NAV
from PPC     import PPC
from SIM_par import SIM_par
from Utils   import spice
from Utils   import events
from Utils   import integrator

# What should run for each simulation execution
def run(SIM_par_override=None, par_override=None, run_id="0000", log_time=True):
    start_time = time.perf_counter()

    # Change default parameters of the simulation
    if SIM_par_override is not None:
        SIM_par.update(SIM_par_override)

    # Load data or run new simulation
    log_save = SIM_par["log_save"]
    log_load = SIM_par["log_load"]
    log_save_path = SIM_par["log_save_path"]
    log_load_path = SIM_par["log_load_path"]

    # Simulation parameters
    sim_dt         = SIM_par["dt"]
    sim_time_start = SIM_par["time_start"]
    sim_time_end   = SIM_par["time_end"]
    sim_time       = sim_time_start
    n_steps        = int((sim_time_end - sim_time_start)/sim_dt) + 1

    # Initialize inputs table
    events_sequence = events.events_sequence
    events_table = events.build_events_table(sim_time_start, sim_time_end, sim_dt, events=events_sequence)

    # Initialize Level-1 modules
    DYN_obj = DYN(         par_override=par_override)
    SEN_obj = SEN(DYN_obj, par_override=par_override)
    NAV_obj = NAV(SEN_obj, par_override=par_override)

    # Initialize timeline
    timeline = {
        "DYN" : PPC.init_timeline(DYN_obj.snapshot(), n_steps),
        "SEN" : PPC.init_timeline(SEN_obj.snapshot(), n_steps),
        "NAV" : PPC.init_timeline(NAV_obj.snapshot(), n_steps),
    }

    # Load timeline from log
    timeline = PPC.load_timeline(timeline, log_load, log_load_path+"/timeline.npz")

    # Main loop
    for step in range(1, n_steps):
        # Update time
        sim_time += sim_dt

        # Load external inputs
        inputs = events_table[sim_time]

        # DYN
        if "DYN" in log_load:
            state = PPC.load_module(timeline["DYN"], step)
            DYN_obj.load_snapshot(state)
        else:
            DYN_obj = integrator.rk4_step(sim_time, sim_dt, DYN_obj, inputs)
            DYN_obj.update_algebraic(sim_time, None, inputs)

        # SEN
        if "SEN" in log_load:
            state = PPC.load_module(timeline["SEN"], step)
            SEN_obj.load_snapshot(state)
        else:
            SEN_obj.update_algebraic(sim_time, DYN_obj, inputs)

        # NAV
        if "NAV" in log_load:
            state = PPC.load_module(timeline["NAV"], step)
            NAV_obj.load_snapshot(state)
        else:
            NAV_obj = integrator.rk4_step(sim_time, sim_dt, NAV_obj, inputs)
            NAV_obj.update_algebraic(sim_time, SEN_obj, inputs)

        # Save results into the timeline
        states = {
            "DYN" : DYN_obj.snapshot(),
            "SEN" : SEN_obj.snapshot(),
            "NAV" : NAV_obj.snapshot(),
        }

        PPC.update_timeline(timeline, states, step)

    # Save log
    if log_save and log_save_path is not None:
        log_save_path = log_save_path+f"_{run_id}"
        PPC.store_timeline(timeline, log_save_path)

    end_time = time.perf_counter()
    if log_time:
        runtime = end_time - start_time
        print(f"Run Time: {runtime:.3f} s")

    return timeline, DYN_obj, SEN_obj, NAV_obj

# What should run after the simulation ends
def stop():
    spice.clear_kernels()
