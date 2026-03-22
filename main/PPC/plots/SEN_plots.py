# List all plots to be generated for SEN modules, organized by group name, to be selected on SIM_par

from PPC import PPC
from Utils import quaternions
from Utils.constants import CONSTANTS_par

def SEN_TIME_plot(timeline, DYN_obj, SEN_obj, NAV_obj):
    time_SIM = timeline["DYN"]["DYN_TIME"]["time_SIM"]
    time_TDB = timeline["DYN"]["DYN_TIME"]["time_TDB"]
    time_OBT = timeline["SEN"]["SEN_TIME"]["time_OBT"]
    fig, ax = PPC.plot(time_SIM, time_TDB, xlabel="time_SIM [s]", ylabel="[s]", label="time_TDB", title="SEN_TIME")
    PPC.plot(time_SIM, time_OBT, label="time_OBT", fig=fig, ax=ax)

def SEN_STR_plot(timeline, DYN_obj, SEN_obj, NAV_obj):
    time_SIM     = timeline["DYN"]["DYN_TIME"]["time_SIM"]
    time_TDB     = timeline["DYN"]["DYN_TIME"]["time_TDB"]
    BOFq_SSB     = timeline["DYN"]["DYN_ATT"]["BOFq_SSB"]
    STARSdir_SSB = timeline["DYN"]["DYN_STR"]["STARSdir_SSB"][0] # [time, star, direction]

    time_STR      = timeline["SEN"]["SEN_STR"]["time_STR"]
    SEN_STRoutflg = timeline["SEN"]["SEN_STR"]["SEN_STRoutflg"]
    field_of_view = SEN_obj.SEN_STR.par["field_of_view"]
    STRq_BOF      = SEN_obj.SEN_STR.par["STRq_BOF"]

    SUNdir_STR_mes    = timeline["SEN"]["SEN_STR"]["SUNdir_STR_mes"]    # [time, direction]
    EARTHdir_STR_mes  = timeline["SEN"]["SEN_STR"]["EARTHdir_STR_mes"]  # [time, direction]
    MOONdir_STR_mes   = timeline["SEN"]["SEN_STR"]["MOONdir_STR_mes"]   # [time, direction]
    MARSdir_STR_mes   = timeline["SEN"]["SEN_STR"]["MARSdir_STR_mes"]   # [time, direction]
    DEIMOSdir_STR_mes = timeline["SEN"]["SEN_STR"]["DEIMOSdir_STR_mes"] # [time, direction]
    PHOBOSdir_STR_mes = timeline["SEN"]["SEN_STR"]["PHOBOSdir_STR_mes"] # [time, direction]
    STARSdir_STR_mes  = timeline["SEN"]["SEN_STR"]["STARSdir_STR_mes"]  # [time, star, direction]

    SUNdir_SC_mes    = timeline["SEN"]["SEN_STR"]["SUNdir_SC_mes"]    # [time, direction]
    EARTHdir_SC_mes  = timeline["SEN"]["SEN_STR"]["EARTHdir_SC_mes"]  # [time, direction]
    MOONdir_SC_mes   = timeline["SEN"]["SEN_STR"]["MOONdir_SC_mes"]   # [time, direction]
    MARSdir_SC_mes   = timeline["SEN"]["SEN_STR"]["MARSdir_SC_mes"]   # [time, direction]
    DEIMOSdir_SC_mes = timeline["SEN"]["SEN_STR"]["DEIMOSdir_SC_mes"] # [time, direction]
    PHOBOSdir_SC_mes = timeline["SEN"]["SEN_STR"]["PHOBOSdir_SC_mes"] # [time, direction]
    STARSdir_SC_mes  = timeline["SEN"]["SEN_STR"]["STARSdir_SC_mes"]  # [time, star, direction]

    # Reshape vectors for plot
    STARSdir_STR_mes = STARSdir_STR_mes.reshape(-1, 3) # [time * star, direction]
    STARSdir_SC_mes  = STARSdir_SC_mes.reshape(-1, 3)  # [time * star, direction]
    STARSdir_SSB     = STARSdir_SSB.reshape(-1, 3)     # [time * star, direction]

    # Plot status
    fig, ax = PPC.plot(time_SIM, SEN_STRoutflg, xlabel="time_SIM [s]", ylabel="SEN_STRoutflg", label="STRoutflag", title="STR output flag")

    # Compare STR and TDB times
    fig, ax = PPC.plot(time_SIM, time_TDB, xlabel="time_SIM [s]", ylabel="time [s]", label="time_TDB", title="STR Time")
    PPC.plot(time_SIM, time_STR, ylabel="time [s]", label="time_STR", fig=fig, ax=ax)

    # Plot STR Gnomonic lens projection
    SUNproj    = PPC.gnomonic_projection(SUNdir_STR_mes)
    EARTHproj  = PPC.gnomonic_projection(EARTHdir_STR_mes)
    MOONproj   = PPC.gnomonic_projection(MOONdir_STR_mes)
    MARSproj   = PPC.gnomonic_projection(MARSdir_STR_mes)
    DEIMOSproj = PPC.gnomonic_projection(DEIMOSdir_STR_mes)
    PHOBOSproj = PPC.gnomonic_projection(PHOBOSdir_STR_mes)
    STARSproj  = PPC.gnomonic_projection(STARSdir_STR_mes)
    boundary   = PPC.gnomonic_boundary(field_of_view)

    fig, ax = PPC.plot([], [], xlabel="X STR", ylabel="Y STR", style='.', label="Stars", aspect="equal", color=PPC.colors["green"])
    PPC.plot(boundary[0], boundary[1], label="FOV", style='--', fig=fig, ax=ax, color=PPC.colors["black"])
    if not PPC.is_nan(SUNproj):    PPC.plot(SUNproj[0],    SUNproj[1],    label="Sun",    style='.',  fig=fig, ax=ax, color=PPC.colors["orange"])
    if not PPC.is_nan(MOONproj):   PPC.plot(MOONproj[0],   MOONproj[1],   label="Moon",   style='.',  fig=fig, ax=ax, color=PPC.colors["grey"])
    if not PPC.is_nan(EARTHproj):  PPC.plot(EARTHproj[0],  EARTHproj[1],  label="Earth",  style='.',  fig=fig, ax=ax, color=PPC.colors["blue"])
    if not PPC.is_nan(DEIMOSproj): PPC.plot(DEIMOSproj[0], DEIMOSproj[1], label="Deimos", style='.',  fig=fig, ax=ax, color=PPC.colors["lightgrey"])
    if not PPC.is_nan(PHOBOSproj): PPC.plot(PHOBOSproj[0], PHOBOSproj[1], label="Phobos", style='.',  fig=fig, ax=ax, color=PPC.colors["darkgrey"])
    if not PPC.is_nan(MARSproj):   PPC.plot(MARSproj[0],   MARSproj[1],   label="Mars",   style='.',  fig=fig, ax=ax, color=PPC.colors["red"])
    if not PPC.is_nan(STARSproj):  PPC.plot(STARSproj[0],  STARSproj[1],                  style='.',  fig=fig, ax=ax, color=PPC.colors["green"])

    # 3D sky sphere
    fig, ax = PPC.plot(STARSdir_SSB[:,0], STARSdir_SSB[:,1], STARSdir_SSB[:,2], style='.', label="Stars", xlabel="X SC", ylabel="Y SC", zlabel="Z SC", title="Star Field Normalized in SC frame", aspect="equal", color=PPC.colors["black"])
    if not PPC.is_nan(STARSdir_SC_mes):   PPC.plot(STARSdir_SC_mes[:,0],  STARSdir_SC_mes[:,1],  STARSdir_SC_mes[:,2],  style='.', label="Visible Stars", fig=fig, ax=ax, color=PPC.colors["green"])
    if not PPC.is_nan(SUNdir_SC_mes):     PPC.plot(SUNdir_SC_mes[:,0],    SUNdir_SC_mes[:,1],    SUNdir_SC_mes[:,1],    style='.', label="Sun",           fig=fig, ax=ax, color=PPC.colors["orange"])
    if not PPC.is_nan(MOONdir_SC_mes):    PPC.plot(MOONdir_SC_mes[:,0],   MOONdir_SC_mes[:,1],   MOONdir_SC_mes[:,2],   style='.', label="MOON",          fig=fig, ax=ax, color=PPC.colors["grey"])
    if not PPC.is_nan(EARTHdir_SC_mes):   PPC.plot(EARTHdir_SC_mes[:,0],  EARTHdir_SC_mes[:,1],  EARTHdir_SC_mes[:,2],  style='.', label="Earth",         fig=fig, ax=ax, color=PPC.colors["blue"])
    if not PPC.is_nan(DEIMOSdir_SC_mes):  PPC.plot(DEIMOSdir_SC_mes[:,0], DEIMOSdir_SC_mes[:,1], DEIMOSdir_SC_mes[:,2], style='.', label="Deimos",        fig=fig, ax=ax, color=PPC.colors["lightgrey"])
    if not PPC.is_nan(PHOBOSdir_SC_mes):  PPC.plot(PHOBOSdir_SC_mes[:,0], PHOBOSdir_SC_mes[:,1], PHOBOSdir_SC_mes[:,2], style='.', label="Phobos",        fig=fig, ax=ax, color=PPC.colors["darkgrey"])
    if not PPC.is_nan(MARSdir_SC_mes):    PPC.plot(MARSdir_SC_mes[:,0],   MARSdir_SC_mes[:,1],   MARSdir_SC_mes[:,2],   style='.', label="Mars",          fig=fig, ax=ax, color=PPC.colors["red"])
    PPC.plot([0], [0], [0],  style='+', label="SC",  fig=fig, ax=ax, color=PPC.colors["magenta"])

    # Get axes in the SSB frame
    STRq_SSB = quaternions.qprod(STRq_BOF, BOFq_SSB)
    BOFx, BOFy, BOFz = quaternions.q2axis(BOFq_SSB)
    STRx, STRy, STRz = quaternions.q2axis(STRq_SSB)

    BOFx = BOFx[-1]
    BOFy = BOFy[-1]
    BOFz = BOFz[-1]

    PPC.plot([0, BOFx[0]], [0, BOFx[1]], [0, BOFx[2]], fig=fig, ax=ax, color=PPC.colors["blue"],  style='--', label="BOFx")
    PPC.plot([0, BOFy[0]], [0, BOFy[1]], [0, BOFy[2]], fig=fig, ax=ax, color=PPC.colors["red"],   style='--', label="BOFy")
    PPC.plot([0, BOFz[0]], [0, BOFz[1]], [0, BOFz[2]], fig=fig, ax=ax, color=PPC.colors["green"], style='--', label="BOFz")

    STRx = STRx[-1]
    STRy = STRy[-1]
    STRz = STRz[-1]

    PPC.plot([0, STRx[0]], [0, STRx[1]], [0, STRx[2]], fig=fig, ax=ax, color=PPC.colors["blue"],  label="STRx")
    PPC.plot([0, STRy[0]], [0, STRy[1]], [0, STRy[2]], fig=fig, ax=ax, color=PPC.colors["red"],   label="STRy")
    PPC.plot([0, STRz[0]], [0, STRz[1]], [0, STRz[2]], fig=fig, ax=ax, color=PPC.colors["green"], label="STRz")

def SEN_PSR_plot(timeline, DYN_obj, SEN_obj, NAV_obj):
    time_SIM  = timeline["DYN"]["DYN_TIME"]["time_SIM"]
    OBTdt_TDB = timeline["DYN"]["DYN_PSR"]["OBTdt_TDB"]

    # TODO: Plot the other SEN_PSR outputs
    time_PSR      = timeline["SEN"]["SEN_PSR"]["time_PSR"]
    SEN_PSRoutflg = timeline["SEN"]["SEN_PSR"]["SEN_PSRoutflg"]
    OBTdt_TDB_mes = timeline["SEN"]["SEN_PSR"]["OBTdt_TDB_mes"]

    # Pulsar signals on the SC
    PULSARname = DYN_obj.DYN_PSR.par["name"]
    n_pulsars  = SEN_obj.SEN_PSR.par["n_pulsars"]

    light_speed_cst = CONSTANTS_par["light_speed_cst"]

    # Measured vs True delay for TOAs between SC and SSB
    fig, ax1 = PPC.plot([], [], ylabel="OBTdt_TDB_mes [s]", title="Measured vs True TOA delay on SC", subplot=(3,1,1))
    fig, ax2 = PPC.plot([], [], ylabel="Time noise [s]", fig=fig, subplot=(3,1,2))
    fig, ax3 = PPC.plot([], [], xlabel="time_SIM [s]", ylabel="Range noise [km]", fig=fig, subplot=(3,1,3))
    time_noise = OBTdt_TDB - OBTdt_TDB_mes
    range_noise = light_speed_cst*time_noise
    for i in range(n_pulsars):
        PPC.plot(time_SIM, OBTdt_TDB_mes[:,i], label=f"{PULSARname[i]} meas.",       fig=fig, ax=ax1)
        PPC.plot(time_SIM, time_noise[:,i],    label=f"{PULSARname[i]} time noise",  fig=fig, ax=ax2)
        PPC.plot(time_SIM, range_noise[:,i],   label=f"{PULSARname[i]} range noise", fig=fig, ax=ax3)

SEN_plots = {
    "SEN_TIME" : SEN_TIME_plot,
    "SEN_STR"  : SEN_STR_plot,
    "SEN_PSR"  : SEN_PSR_plot,
}
