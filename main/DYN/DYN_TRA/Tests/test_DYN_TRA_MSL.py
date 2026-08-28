# DYN_TRA Module Test
# Use the main code solver to propagate spacecraft trajectory
# Compares the result using data from Spice
# Run from repository root with:
#   python3 -m DYN.DYN_TRA.Tests.test_DYN_TRA_MSL

import numpy as np
from DYN.DYN import DYN
from PPC     import PPC
from Utils   import spice
from Utils   import events
from Utils   import integrator

##################
# MSL PARAMETERS #
##################

# Load NASA's MSL data with SPICE
# https://naif.jpl.nasa.gov/pub/naif/MSL/kernels/
kernel_dir = "Utils/kernels/"
spice.load_kernel(kernel_dir + "msl_cruise.bsp")

# Get Spacecraft initial condition
# MSL trajectory will be compared during the cruise phase, in a period of time without any trajectory correction maneuvers
time_TDB_ini = spice.get_time("2012-01-01 T00:00:00")
time_TDB_end = spice.get_time("2012-07-30 T00:00:00")

MSLpos_SSB_ini, MSLvel_SSB_ini = spice.get_state("MSL", time_TDB_ini)

#########
# SETUP #
#########

SIM_par = {
    # Time parameters
    "dt"           : 600, # [s] 10 min
    "time_start"   : 0,   # [s]
    "time_end"     : int(time_TDB_end - time_TDB_ini), # [s]

    # Log parameters
    "DYN_log_save" : False,
    "DYN_log_load" : False,
    "DYN_log_path" : "Logs/DYN",

    # List of default plots
    "PPC_plot_list" : [
        "DYN_TIME",
        "DYN_EPH",
        "DYN_TRA",
    ],
}

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

######################
# INITIALIZE MODULES #
######################

DYN_TIME_par = {
    "time_SIM_ini" : 0,
    "time_TDB_ini"  : time_TDB_ini,
}

DYN_TRA_par = {
    "ref_elements"  : "rvi",
    "BODY_ini"      : "SUN",
    "sma_ini"       : None,
    "ecc_ini"       : None,
    "incl_ini"      : None,
    "raan_ini"      : None,
    "argp_ini"      : None,
    "tano_ini"      : None,
    "SCpos_SSB_ini" : MSLpos_SSB_ini, # [km]   Initial position relative to SSB frame
    "SCvel_SSB_ini" : MSLvel_SSB_ini, # [km/s] Initial velocity relative to SSB frame
}

par_override = {
    "DYN_TIME" : DYN_TIME_par,
    "DYN_TRA"  : DYN_TRA_par,
}

# Normally initialize DYN module, with override parameters
DYN_obj = DYN(par_override)
# Set only modules needed for the test
DYN_obj.modules = [
    DYN_obj.DYN_TIME,
    DYN_obj.DYN_EPH,
    DYN_obj.DYN_EPH,
    DYN_obj.DYN_EPH,
    DYN_obj.DYN_TRA,
]

# Define MSL structure to store its states
MSL_TRA = {
    "MSLpos_SSB" : MSLpos_SSB_ini,
    "MSLvel_SSB" : MSLvel_SSB_ini,
}

# Initialize timeline
timeline = {"DYN" : PPC.init_timeline(DYN_obj.snapshot(), n_steps),
            "MSL" : PPC.init_timeline(MSL_TRA, n_steps)}

########
# LOOP #
########

# Main loop
for step in range(1, n_steps):
    # Update time
    sim_time += sim_dt

    # Load external inputs
    inputs = events_table[sim_time]

    # Integrate all dynamic states together
    DYN_obj = integrator.rk4_step(sim_time, sim_dt, DYN_obj, inputs)
    # Update algebraic modules
    DYN_obj.update_algebraic(sim_time, None, inputs)

    # Update MSL state
    time_TDB = DYN_obj.snapshot()["DYN_TIME"]["time_TDB"]
    MSLpos_SSB, MSLvel_SSB = spice.get_state("MSL", time_TDB)

    # Save results into the timeline
    output = {"DYN" : DYN_obj.snapshot(),
                "MSL" : {"MSLpos_SSB" : MSLpos_SSB, "MSLvel_SSB" : MSLvel_SSB}}
    PPC.update_timeline(timeline, output, step)

############
# TEARDOWN #
############

# Clear Kernels from memory
spice.clear_kernels()

###########
# RESULTS #
###########

SCpos_SSB    = timeline["DYN"]["DYN_TRA"]["SCpos_SSB"]
SCvel_SSB    = timeline["DYN"]["DYN_TRA"]["SCvel_SSB"]
MSLpos_SSB   = timeline["MSL"]["MSLpos_SSB"]
MSLvel_SSB   = timeline["MSL"]["MSLvel_SSB"]
EARTHpos_SSB = timeline["DYN"]["DYN_EPH"]["EARTHpos_SSB"]
MARSpos_SSB  = timeline["DYN"]["DYN_EPH"]["MARSpos_SSB"]
MARSvel_SSB  = timeline["DYN"]["DYN_EPH"]["MARSvel_SSB"]
SUNpos_SSB   = timeline["DYN"]["DYN_EPH"]["SUNpos_SSB"]

time_days = timeline["DYN"]["DYN_TIME"]["time_SIM"]/3600/24

# Plot spacecraft position vs MSL position
fig, ax_x = PPC.plot(time_days, SCpos_SSB[:,0], ylabel="posx_SSB [km]", label="SC", title="pos SC vs MSL", subplot=(3,1,1))
fig, ax_y = PPC.plot(time_days, SCpos_SSB[:,1], ylabel="posy_SSB [km]", label="SC", fig=fig, subplot=(3,1,2))
fig, ax_z = PPC.plot(time_days, SCpos_SSB[:,2], xlabel="Time SIM [days]", ylabel="posz_SSB [km]", label="SC", fig=fig, subplot=(3,1,3))
PPC.plot(time_days, MSLpos_SSB[:,0], label="MSL", fig=fig, ax=ax_x)
PPC.plot(time_days, MSLpos_SSB[:,1], label="MSL", fig=fig, ax=ax_y)
PPC.plot(time_days, MSLpos_SSB[:,2], label="MSL", fig=fig, ax=ax_z)

# Plot spacecraft velocity vs MSL velocity
fig, ax_x = PPC.plot(time_days, SCvel_SSB[:,0], ylabel="velx_SSB [km/s]", label="SC", title="vel SC vs MSL", subplot=(3,1,1))
fig, ax_y = PPC.plot(time_days, SCvel_SSB[:,1], ylabel="vely_SSB [km/s]", label="SC", fig=fig, subplot=(3,1,2))
fig, ax_z = PPC.plot(time_days, SCvel_SSB[:,2], xlabel="Time SIM [days]", ylabel="velz_SSB [km/s]", label="SC", fig=fig, subplot=(3,1,3))
PPC.plot(time_days, MSLvel_SSB[:,0], label="MSL", fig=fig, ax=ax_x)
PPC.plot(time_days, MSLvel_SSB[:,1], label="MSL", fig=fig, ax=ax_y)
PPC.plot(time_days, MSLvel_SSB[:,2], label="MSL", fig=fig, ax=ax_z)

# Plot error on position comparing SC and MSL
pos_error_x = np.sqrt((SCpos_SSB[:,0] - MSLpos_SSB[:,0])**2)
pos_error_y = np.sqrt((SCpos_SSB[:,1] - MSLpos_SSB[:,1])**2)
pos_error_z = np.sqrt((SCpos_SSB[:,2] - MSLpos_SSB[:,2])**2)
pos_error = np.sqrt(pos_error_x**2 + pos_error_y**2 + pos_error_z**2)
fig, ax_r = PPC.plot(time_days, pos_error_x, ylabel="pos error [km]", label="|$r_{SC,x} - r_{MSL,x}$|", title="Error pos SC vs MSL", subplot=(2,1,1), color=PPC.colors['blue'])
PPC.plot(time_days, pos_error_y, label="|$r_{SC,y} - r_{MSL,y}$|", fig=fig, ax=ax_r, color=PPC.colors['red'])
PPC.plot(time_days, pos_error_z, label="|$r_{SC,z} - r_{MSL,z}$|", fig=fig, ax=ax_r, color=PPC.colors['green'])
PPC.plot(time_days, pos_error, label="|$r_{SC} - r_{MSL}$|", fig=fig, ax=ax_r, color=PPC.colors['magenta'])

# Plot error on velocity comparing SC and MSL
vel_error_x = np.sqrt((SCvel_SSB[:,0] - MSLvel_SSB[:,0])**2)
vel_error_y = np.sqrt((SCvel_SSB[:,1] - MSLvel_SSB[:,1])**2)
vel_error_z = np.sqrt((SCvel_SSB[:,2] - MSLvel_SSB[:,2])**2)
vel_error = np.sqrt(vel_error_x**2 + vel_error_y**2 + vel_error_z**2)
fig, ax_v = PPC.plot(time_days, vel_error_x, ylabel="vel error [km/s]", label="|$v_{SC,x} - v_{MSL,x}$|", title="Error vel SC vs MSL", subplot=(2,1,2), color=PPC.colors['blue'], fig=fig)
PPC.plot(time_days, vel_error_y, label="|$v_{SC,y} - v_{MSL,y}$|", subplot=(4,1,2), fig=fig, ax=ax_v, color=PPC.colors['red'])
PPC.plot(time_days, vel_error_z, label="|$v_{SC,z} - v_{MSL,z}$|", subplot=(4,1,3), fig=fig, ax=ax_v, color=PPC.colors['green'])
PPC.plot(time_days, vel_error, xlabel="Time SIM [days]", ylabel="vel error [km/s]", label="|$v_{SC} - v_{MSL}$|", fig=fig, ax=ax_v, color=PPC.colors['magenta'])

# Trajectoy in the orbital plane
fig, ax = PPC.plot(SCpos_SSB[:,0], SCpos_SSB[:,1], xlabel="x SSB [km]", ylabel="y SSB [km]", title="SC vs MSL trajectory", aspect='equal', color=PPC.colors["magenta"], style=':', zorder=4)
PPC.plot(SCpos_SSB[-1,0],    SCpos_SSB[-1,1],    label='SC', style='x', fig=fig, ax=ax, color=PPC.colors["magenta"], zorder=4)
PPC.plot(MSLpos_SSB[:,0],    MSLpos_SSB[:,1],    fig=fig, ax=ax, color=PPC.colors["grey"])
PPC.plot(MSLpos_SSB[-1,0],   MSLpos_SSB[-1,1],   label="MSL", style='X', fig=fig, ax=ax, color=PPC.colors["grey"])
PPC.plot(EARTHpos_SSB[:,0],  EARTHpos_SSB[:,1],  fig=fig, ax=ax, color=PPC.colors["blue"])
PPC.plot(EARTHpos_SSB[-1,0], EARTHpos_SSB[-1,1], label="Earth", style='o', fig=fig, ax=ax, color=PPC.colors["blue"])
PPC.plot(MARSpos_SSB[:,0],   MARSpos_SSB[:,1],   fig=fig, ax=ax, color=PPC.colors["red"])
PPC.plot(MARSpos_SSB[-1,0],  MARSpos_SSB[-1,1],  label="Mars", style='o', fig=fig, ax=ax, color=PPC.colors["red"])
PPC.plot(SUNpos_SSB[0,0],    SUNpos_SSB[0,1],    label="Sun", style='o', fig=fig, ax=ax, color=PPC.colors["orange"])

# 3D trajectories
fig, ax = PPC.plot([0], [0], [0], style='+', label="SSB", xlabel="X SSB [km]", ylabel="Y SSB [km]", zlabel="Z SSB [km]", title="Simulation Trajectories", aspect="equal",                      color=PPC.colors["black"])
bodies = [
    dict(name="SUN",    idx=0, color=PPC.colors["orange"]),
    dict(name="MOON",   idx=2, color=PPC.colors["grey"]),
    dict(name="EARTH",  idx=1, color=PPC.colors["blue"]),
    dict(name="DEIMOS", idx=4, color=PPC.colors["darkgrey"]),
    dict(name="PHOBOS", idx=5, color=PPC.colors["lightgrey"]),
    dict(name="MARS",   idx=3, color=PPC.colors["red"]),
]
for body in bodies:
    name = body["name"]
    PPC.plot(timeline["DYN"]["DYN_EPH"][f"{name}pos_SSB"][:,0], timeline["DYN"]["DYN_EPH"][f"{name}pos_SSB"][:,1], timeline["DYN"]["DYN_EPH"][f"{name}pos_SSB"][:,2], label=name, fig=fig, ax=ax, color=body["color"])
# Plot Spacecraft
PPC.plot(SCpos_SSB[:,0],  SCpos_SSB[:,1],  SCpos_SSB[:,2],  label="SC",  fig=fig, ax=ax, color=PPC.colors["magenta"])
PPC.plot(MSLpos_SSB[:,0], MSLpos_SSB[:,1], MSLpos_SSB[:,2], label="MSL", fig=fig, ax=ax, color=PPC.colors["grey"])

# SC and MSL distance to Mars
SCpos_MCI  = SCpos_SSB  - MARSpos_SSB
MSLpos_MCI = MSLpos_SSB - MARSpos_SSB
SCvel_MCI  = SCvel_SSB  - MARSvel_SSB
MSLvel_MCI = MSLvel_SSB - MARSvel_SSB

SCpos_error  = np.sqrt(SCpos_MCI[:,0]**2 + SCpos_MCI[:,1]**2 + SCpos_MCI[:,2]**2)
MSLpos_error = np.sqrt(MSLpos_MCI[:,0]**2 + MSLpos_MCI[:,1]**2 + MSLpos_MCI[:,2]**2)
fig, ax_r = PPC.plot(time_days, SCpos_error, ylabel="SC distance to Mars [km]", label="|$r_{SC}^{MCI}$|", title="SC distance to Mars", subplot=(2,1,1), color=PPC.colors['black'])
PPC.plot(time_days, MSLpos_error, label="|$r_{MSL}^{MCI}$|", fig=fig, ax=ax_r, color=PPC.colors['green'])

SCvel_error  = np.sqrt(SCvel_MCI[:,0]**2 + SCvel_MCI[:,1]**2 + SCvel_MCI[:,2]**2)
MSLvel_error = np.sqrt(MSLvel_MCI[:,0]**2 + MSLvel_MCI[:,1]**2 + MSLvel_MCI[:,2]**2)
fig, ax_v = PPC.plot(time_days, SCvel_error, ylabel="SC velocity to Mars [km/s]", label="|$v_{SC}^{MCI}$|", subplot=(2,1,2), color=PPC.colors['black'], fig=fig)
PPC.plot(time_days, MSLvel_error, label="|$v_{MSL}^{MCI}$|", fig=fig, ax=ax_v, color=PPC.colors['green'])

PPC.show_plot()

# Save data
save_flg = True
if save_flg:
    # Downsample data
    n_samples = 200
    SCpos_SSB_ds    = PPC.downsample(SCpos_SSB,    n_samples)
    SCvel_SSB_ds    = PPC.downsample(SCvel_SSB,    n_samples)
    MSLpos_SSB_ds   = PPC.downsample(MSLpos_SSB,   n_samples)
    MSLvel_SSB_ds   = PPC.downsample(MSLvel_SSB,   n_samples)
    EARTHpos_SSB_ds = PPC.downsample(EARTHpos_SSB, n_samples)
    MARSpos_SSB_ds  = PPC.downsample(MARSpos_SSB,  n_samples)
    SUNpos_SSB_ds   = PPC.downsample(SUNpos_SSB,   n_samples)
    pos_error_ds    = PPC.downsample(pos_error,    n_samples)
    vel_error_ds    = PPC.downsample(vel_error,    n_samples)

    # Write data into file
    np.savetxt("SCpos_SSB.csv",    SCpos_SSB_ds,    delimiter=',')
    np.savetxt("SCvel_SSB.csv",    SCvel_SSB_ds,    delimiter=',')
    np.savetxt("MSLpos_SSB.csv",   MSLpos_SSB_ds,   delimiter=',')
    np.savetxt("MSLvel_SSB.csv",   MSLvel_SSB_ds,   delimiter=',')
    np.savetxt("EARTHpos_SSB.csv", EARTHpos_SSB_ds, delimiter=',')
    np.savetxt("MARSpos_SSB.csv",  MARSpos_SSB_ds,  delimiter=',')
    np.savetxt("SUNpos_SSB.csv",   SUNpos_SSB_ds,   delimiter=',')
    np.savetxt("pos_error.csv",    pos_error_ds,    delimiter=',')
    np.savetxt("vel_error.csv",    vel_error_ds,    delimiter=',')
