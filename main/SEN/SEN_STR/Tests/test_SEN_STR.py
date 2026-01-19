# SEN_STR Module Test
# Generates the sensor field of view when the spacecraft is aways pointing to Mars
# Compares the measurements with the simulated stars from DYN_STR
# Run from repository root with:
#   python3 -m SEN.SEN_STR.Tests.test_SEN_STR

import numpy as np
from DYN.DYN       import DYN
from SEN.SEN       import SEN
from PPC           import PPC
from PPC.PPC_plots import PPC_plots
from PPC.PPC_plots import colors
from Utils         import spice
from Utils         import events
from Utils         import integrator
from Utils         import quaternions

#########
# SETUP #
#########

SIM_par = {
    # Time parameters
    "dt"           : 600,       # [s] 10 min
    "time_start"   : 0,         # [s]
    "time_end"     : 3600*24*5, # [s] 5 days

    # Log parameters
    "DYN_log_save" : False,
    "DYN_log_load" : False,
    "DYN_log_path" : "Logs/DYN",

    # List of default plots
    "PPC_plot_list" : [
        # "SEN_STR",
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

# Set events
events_sequence = events.events_sequence
events_sequence["SEN.SEN_STR.STRenableflg"] = [
        (100,   0),
        (6000,  1),
        (12000, 0),
        (13000, 1),
]

# Initialize inputs table
events_table = events.build_events_table(sim_time_start, sim_time_end, sim_dt, events=events_sequence)

######################
# INITIALIZE MODULES #
######################

# Normally initialize DYN modules
DYN_obj = DYN()

# Set only modules needed for the test
DYN_obj.modules = [
    DYN_obj.DYN_TIME,
    DYN_obj.DYN_SUN,
    DYN_obj.DYN_EARTH,
    DYN_obj.DYN_MARS,
    DYN_obj.DYN_TRA,
    DYN_obj.DYN_GRV,
    DYN_obj.DYN_ATT,
    DYN_obj.DYN_STR,
]

# Normally initialize SEN modules
SEN_obj = SEN(DYN_obj)

# Set only modules needed for the test
SEN_obj.modules = [
    SEN_obj.SEN_STR,
]

# Initialize timeline
timeline = {
    "DYN" : PPC.init_timeline(DYN_obj.snapshot(), n_steps),
    "SEN" : PPC.init_timeline(SEN_obj.snapshot(), n_steps)
}

########
# LOOP #
########

# Main loop
for step in range(1, n_steps):

    # Load external inputs
    inputs = events_table[sim_time]

    ###################
    # UPDATE ATTITUDE #
    ###################

    # Spacecraft is aways pointing +STRz to Mars
    DYN_last_state = DYN_obj.snapshot()
    MARSpos_SSB = DYN_last_state["DYN_MARS"]["MARSpos_SSB"]
    SCpos_SSB   = DYN_last_state["DYN_TRA"]["SCpos_SSB"]
    STRy_STR    = np.array([0,1,0])
    STRz_STR    = np.array([0,0,1])
    STRq_BOF    = SEN_obj.SEN_STR.par["STRq_BOF"]

    BOFq_STR    = quaternions.qtrans(STRq_BOF)

    # Line of sight from SC to Mars
    MARSpos_SC = MARSpos_SSB - SCpos_SSB
    MARSdir_SC = SEN_obj.SEN_STR.dir_from_pos(MARSpos_SC)

    # STR boresight in the BOF frame
    STRz_BOF = quaternions.qvecprod(BOFq_STR, STRz_STR)

    # Get the quaternion which points STRz_BOF to MARSu_SSB
    BOFq_SSB_tgt = quaternions.vecqvec(STRz_BOF, MARSdir_SC)

    inputs["DYN"]["DYN_ATT"]["BOFq_SSB"] = BOFq_SSB_tgt

    ###################

    # Update algebraic modules first
    DYN_obj.update_algebraic(sim_time, inputs)
    # Integrate all dynamic states together
    integrator.rk4_step(sim_time, sim_dt, DYN_obj, inputs)
    # Update SEN
    SEN_obj.update_algebraic(sim_time, DYN_obj, inputs)

    # Save results into the timeline
    output = {
        "DYN" : DYN_obj.snapshot(),
        "SEN" : SEN_obj.snapshot(),
    }
    PPC.update_timeline(timeline, output, step)

    # Update time
    sim_time += sim_dt

############
# TEARDOWN #
############

# Clear Kernels from memory
spice.clear_kernels()

###########
# RESULTS #
###########

field_of_view = SEN_obj.SEN_STR.par["field_of_view"]
STRq_BOF      = SEN_obj.SEN_STR.par["STRq_BOF"]

time_SIM          = timeline["DYN"]["DYN_TIME"]["time_SIM"]
time_STR          = timeline["SEN"]["SEN_STR"]["time_STR"]
STRoutflg         = timeline["SEN"]["SEN_STR"]["STRoutflg"]
SUNdir_STR_mes    = timeline["SEN"]["SEN_STR"]["SUNdir_STR_mes"].T
EARTHdir_STR_mes  = timeline["SEN"]["SEN_STR"]["EARTHdir_STR_mes"].T
MARSdir_STR_mes   = timeline["SEN"]["SEN_STR"]["MARSdir_STR_mes"].T
STARSdir_STR_mes  = timeline["SEN"]["SEN_STR"]["STARSdir_STR_mes"].T # [coord, star, time]
SUNdir_SSB_mes    = timeline["SEN"]["SEN_STR"]["SUNdir_SSB_mes"].T
EARTHdir_SSB_mes  = timeline["SEN"]["SEN_STR"]["EARTHdir_SSB_mes"].T
MARSdir_SSB_mes   = timeline["SEN"]["SEN_STR"]["MARSdir_SSB_mes"].T
STARSdir_SSB_mes  = timeline["SEN"]["SEN_STR"]["STARSdir_SSB_mes"].T
BOFq_SSB_mes      = timeline["SEN"]["SEN_STR"]["BOFq_SSB_mes"]

SCpos_SSB    = timeline["DYN"]["DYN_TRA"]["SCpos_SSB"].T
STARSdir_SSB = timeline["DYN"]["DYN_STR"]["STARSdir_SSB"].T
SUNpos_SSB   = timeline["DYN"]["DYN_SUN"]["SUNpos_SSB"].T
EARTHpos_SSB = timeline["DYN"]["DYN_EARTH"]["EARTHpos_SSB"].T
MARSpos_SSB  = timeline["DYN"]["DYN_MARS"]["MARSpos_SSB"].T

BOFq_SSB = timeline["DYN"]["DYN_ATT"]["BOFq_SSB"]

# Get object directions
SUNpos_SC   = SUNpos_SSB   - SCpos_SSB
EARTHpos_SC = EARTHpos_SSB - SCpos_SSB
MARSpos_SC  = MARSpos_SSB  - SCpos_SSB

SUNdir_SC   = SEN_obj.SEN_STR.dir_from_pos(SUNpos_SC.T).T
EARTHdir_SC = SEN_obj.SEN_STR.dir_from_pos(EARTHpos_SC.T).T
MARSdir_SC  = SEN_obj.SEN_STR.dir_from_pos(MARSpos_SC.T).T

# Plot status
fig, ax = PPC.plot(time_SIM, STRoutflg, ylabel="STRoutflg", label="STRoutflag", title="STR output flag", subplot=(2,1,1))

# Compare STR and SIM times
fig, ax = PPC.plot(time_SIM, time_SIM, xlabel="time_SIM [s]", ylabel="time [s]", label="time_SIM", title="STR Time", fig=fig, subplot=(2,1,2))
PPC.plot(time_SIM, time_STR, ylabel="time [s]", label="time_STR", fig=fig, ax=ax)

# Plot BOFq_SSB
fig, ax0 = PPC.plot(time_SIM, BOFq_SSB[:,0], label="BOFq_SSB[0]", title="Spacecraft orientation in the SSB frame", subplot=(4,1,1))
fig, ax1 = PPC.plot(time_SIM, BOFq_SSB[:,1], label="BOFq_SSB[1]", subplot=(4,1,2), fig=fig)
fig, ax2 = PPC.plot(time_SIM, BOFq_SSB[:,2], label="BOFq_SSB[2]", subplot=(4,1,3), fig=fig)
fig, ax3 = PPC.plot(time_SIM, BOFq_SSB[:,3], xlabel="time_SIM [s]", label="BOFq_SSB[3]", subplot=(4,1,4), fig=fig)

fig, ax0 = PPC.plot(time_SIM, BOFq_SSB_mes[:,0], label="BOFq_SSB_mes[0]", subplot=(4,1,1), fig=fig, ax=ax0)
fig, ax1 = PPC.plot(time_SIM, BOFq_SSB_mes[:,1], label="BOFq_SSB_mes[1]", subplot=(4,1,2), fig=fig, ax=ax1)
fig, ax2 = PPC.plot(time_SIM, BOFq_SSB_mes[:,2], label="BOFq_SSB_mes[2]", subplot=(4,1,3), fig=fig, ax=ax2)
fig, ax3 = PPC.plot(time_SIM, BOFq_SSB_mes[:,3], label="BOFq_SSB_mes[3]", subplot=(4,1,4), fig=fig, ax=ax3)

# Plot STR gnomonic lens projection
SUNproj   = PPC.gnomonic_projection(SUNdir_STR_mes)
EARTHproj = PPC.gnomonic_projection(EARTHdir_STR_mes)
MARSproj  = PPC.gnomonic_projection(MARSdir_STR_mes)
STARSproj = PPC.gnomonic_projection(STARSdir_STR_mes)
boundary  = PPC.gnomonic_boundary(field_of_view)

fig, ax = PPC.plot([], [], xlabel="X STR", ylabel="Y STR", style='.', label="Stars", aspect="equal", color=colors["green"])
# PPC.plot(SUNproj[0],   SUNproj[1],   label="Sun",   style='.',  fig=fig, ax=ax, color=colors["orange"])
# PPC.plot(EARTHproj[0], EARTHproj[1], label="Earth", style='.',  fig=fig, ax=ax, color=colors["blue"])
PPC.plot(MARSproj[0],  MARSproj[1],  label="Mars",  style='.',  fig=fig, ax=ax, color=colors["red"])
PPC.plot(STARSproj[0], STARSproj[1],                style='.',  fig=fig, ax=ax, color=colors["green"])
PPC.plot(boundary[0], boundary[1],   label="FOV",   style='--', fig=fig, ax=ax, color=colors["darkgrey"])

# 3D sky sphere on step 0
fig, ax = PPC.plot(STARSdir_SSB[0,:,-1], STARSdir_SSB[1,:,-1], STARSdir_SSB[2,:,-1], style='.', label="Stars", xlabel="X SSB", ylabel="Y SSB", zlabel="Z SSB", title="Star Field Normalized in SSB frame", aspect="equal", color=colors["black"])
PPC.plot(SUNdir_SC[0,:],     SUNdir_SC[1,:],     SUNdir_SC[2,:],     style='.', label="Sun",           fig=fig, ax=ax, color=colors["orange"])
PPC.plot(EARTHdir_SC[0,:],   EARTHdir_SC[1,:],   EARTHdir_SC[2,:],   style='.', label="Earth",         fig=fig, ax=ax, color=colors["blue"])
PPC.plot(MARSdir_SC[0,:],    MARSdir_SC[1,:],    MARSdir_SC[2,:],    style='.', label="Mars",          fig=fig, ax=ax, color=colors["red"])

PPC.plot(STARSdir_SSB_mes[0,:,-1], STARSdir_SSB_mes[1,:,-1], STARSdir_SSB_mes[2,:,-1], style='x', label="Visible Stars", fig=fig, ax=ax, color=colors["green"])
# PPC.plot(SUNdir_SSB_mes[0,:],     SUNdir_SSB_mes[1,:],     SUNdir_SSB_mes[2,:],     style='x', label="Sun Meas",           fig=fig, ax=ax, color=colors["orange"])
# PPC.plot(EARTHdir_SSB_mes[0,:],   EARTHdir_SSB_mes[1,:],   EARTHdir_SSB_mes[2,:],   style='x', label="Earth Meas",         fig=fig, ax=ax, color=colors["blue"])
PPC.plot(MARSdir_SSB_mes[0,:],    MARSdir_SSB_mes[1,:],    MARSdir_SSB_mes[2,:],    style='x', label="Visible Mars",          fig=fig, ax=ax, color=colors["red"])
PPC.plot([0], [0], [0],  style='+', label="SC",  fig=fig, ax=ax, color=colors["magenta"])

# Get axes in the SSB frame

STRq_SSB = quaternions.qprod(STRq_BOF, BOFq_SSB)
BOFx, BOFy, BOFz = quaternions.q2axis(BOFq_SSB)
STRx, STRy, STRz = quaternions.q2axis(STRq_SSB)

BOFx = BOFx[-1]
BOFy = BOFy[-1]
BOFz = BOFz[-1]

PPC.plot([0, BOFx[0]], [0, BOFx[1]], [0, BOFx[2]], fig=fig, ax=ax, color=colors["blue"],  style='--', label="BOFx")
PPC.plot([0, BOFy[0]], [0, BOFy[1]], [0, BOFy[2]], fig=fig, ax=ax, color=colors["red"],   style='--', label="BOFy")
PPC.plot([0, BOFz[0]], [0, BOFz[1]], [0, BOFz[2]], fig=fig, ax=ax, color=colors["green"], style='--', label="BOFz")

STRx = STRx[-1]
STRy = STRy[-1]
STRz = STRz[-1]

PPC.plot([0, STRx[0]], [0, STRx[1]], [0, STRx[2]], fig=fig, ax=ax, color=colors["blue"],  label="STRx")
PPC.plot([0, STRy[0]], [0, STRy[1]], [0, STRy[2]], fig=fig, ax=ax, color=colors["red"],   label="STRy")
PPC.plot([0, STRz[0]], [0, STRz[1]], [0, STRz[2]], fig=fig, ax=ax, color=colors["green"], label="STRz")

# Plot true vs measured Mars direction over time
fig, ax = PPC.plot(time_SIM, MARSdir_SC.T, label=["x","y","z"], xlabel="time_SIM [s]", ylabel="MARSdir_SSB", title="Mars direction measurement")
PPC.plot(time_SIM, MARSdir_SSB_mes.T, label=["x_mes","y_mes","z_mes"], xlabel="time_SIM [s]", ylabel="MARSdir_SSB_mes", style='--', fig=fig, ax=ax)

for key in SIM_par["PPC_plot_list"]:
    if key in PPC_plots:
        PPC_plots[key](timeline, DYN_obj, SEN_obj)
PPC.show_plot()
