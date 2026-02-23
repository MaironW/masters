# List all plots to be generated, organized by group name, to be selected on SIM_par

import numpy as np
from PPC import PPC
from Utils import quaternions
from Utils.constants import CONSTANTS_par

# Define color cycle
colors = {
    "blue"      : "#0066FF",
    "red"       : "#CC0000",
    "green"     : "#33CC00",
    "magenta"   : "#FF00FF",
    "orange"    : "#FE9920",
    "purple"    : "#8052CF",
    "yellow"    : "#F5F22B",
    "grey"      : "#505050",
    "lightgrey" : "#707070",
    "darkgrey"  : "#303030",
    "black"     : "#111111",
}
PPC.setup_plot(colors.values())

def DYN_TIME_plot(timeline, DYN_obj, SEN_obj, NAV_obj):
    PPC.plot(timeline["DYN"]["DYN_TIME"]["time_SIM"], timeline["DYN"]["DYN_TIME"]["time_TDB"], xlabel="time_SIM [s]", ylabel="time_TDB [s]", title="Time")

def DYN_EPH_plot(timeline, DYN_obj, SEN_obj, NAV_obj):
    time_SIM = timeline["DYN"]["DYN_TIME"]["time_SIM"]

    bodies = [
        dict(name="SUN",    idx=0, color=colors["orange"]),
        dict(name="MOON",   idx=2, color=colors["grey"]),
        dict(name="EARTH",  idx=1, color=colors["blue"]),
        dict(name="DEIMOS", idx=4, color=colors["darkgrey"]),
        dict(name="PHOBOS", idx=5, color=colors["lightgrey"]),
        dict(name="MARS",   idx=3, color=colors["red"]),
    ]

    for body in bodies:
        name = body["name"]
        # Plot reference position from DYN
        fig, ax1 = PPC.plot(time_SIM, timeline["DYN"]["DYN_EPH"][f"{name}pos_SSB"], ylabel=f"{name}pos_SSB [km]", label=["x","y","z"], title=f"{name}pos_SSB", subplot=(2,1,1))
        # Plot reference velocity from DYN
        fig, ax2 = PPC.plot(time_SIM, timeline["DYN"]["DYN_EPH"][f"{name}vel_SSB"], ylabel=f"{name}vel_SSB [km]", label=["x","y","z"], title=f"{name}vel_SSB", fig=fig, subplot=(2,1,2))

def DYN_GRV_plot(timeline, DYN_obj, SEN_obj, NAV_obj):
    # Gravity acceleration on SSB frame
    PPC.plot(timeline["DYN"]["DYN_TIME"]["time_SIM"], timeline["DYN"]["DYN_GRV"]["grvacc_SSB"], ylabel="grvacc_SSB [km/s^2]", label=["x","y","z"], title="grvacc_SSB")

def DYN_TRA_plot(timeline, DYN_obj, SEN_obj, NAV_obj):
    # 3D trajectories
    fig, ax = PPC.plot([0], [0], [0], style='+', label="SSB", xlabel="X SSB [km]", ylabel="Y SSB [km]", zlabel="Z SSB [km]", title="Simulation Trajectories", aspect="equal",                      color=colors["black"])
    bodies = [
        dict(name="SUN",    idx=0, color=colors["orange"]),
        dict(name="MOON",   idx=2, color=colors["grey"]),
        dict(name="EARTH",  idx=1, color=colors["blue"]),
        dict(name="DEIMOS", idx=4, color=colors["darkgrey"]),
        dict(name="PHOBOS", idx=5, color=colors["lightgrey"]),
        dict(name="MARS",   idx=3, color=colors["red"]),
    ]
    for body in bodies:
        name = body["name"]
        PPC.plot(timeline["DYN"]["DYN_EPH"][f"{name}pos_SSB"][:,0], timeline["DYN"]["DYN_EPH"][f"{name}pos_SSB"][:,1], timeline["DYN"]["DYN_EPH"][f"{name}pos_SSB"][:,2], label=name, fig=fig, ax=ax, color=body["color"])
    # Plot Spacecraft
    PPC.plot(timeline["DYN"]["DYN_TRA"]["SCpos_SSB"][:,0], timeline["DYN"]["DYN_TRA"]["SCpos_SSB"][:,1], timeline["DYN"]["DYN_TRA"]["SCpos_SSB"][:,2], label="SC", fig=fig, ax=ax, color=colors["magenta"])

    # 2D SC and Moon trajectories around Earth
    MOONpos_ECI = timeline["DYN"]["DYN_EPH"]["EARTHpos_SSB"] - timeline["DYN"]["DYN_EPH"]["MOONpos_SSB"]
    fig, ax = PPC.plot([0], [0], style='o', label="Earth", xlabel="X ECI [km]", ylabel="Y ECI [km]", title="Trajectories around Earth", aspect="equal", color=colors["blue"])
    PPC.plot(timeline["DYN"]["DYN_TRA"]["SCpos_ECI"][:,0], timeline["DYN"]["DYN_TRA"]["SCpos_ECI"][:,1], label="SC",   fig=fig, ax=ax, color=colors["magenta"])
    PPC.plot(MOONpos_ECI[:,0],                             MOONpos_ECI[:,1],                             label="Moon", fig=fig, ax=ax, color=colors["grey"])

    # Spacecraft state
    fig, ax = PPC.plot(timeline["DYN"]["DYN_TIME"]["time_SIM"], timeline["DYN"]["DYN_TRA"]["SCpos_SSB"], ylabel="SCpos_SSB [km]", label=["x","y","z"], title="SCpos_SSB", subplot=(2,1,1))
    PPC.plot(timeline["DYN"]["DYN_TIME"]["time_SIM"], timeline["DYN"]["DYN_TRA"]["SCvel_SSB"], ylabel="SCvel_SSB [km/s]", label=["x","y","z"], title="SCvel_SSB", fig=fig, subplot=(2,1,2))

def DYN_ATT_plot(timeline, DYN_obj, SEN_obj, NAV_obj):
    # Spacecraft attitude
    PPC.plot(timeline["DYN"]["DYN_TIME"]["time_SIM"], timeline["DYN"]["DYN_ATT"]["BOFq_SSB"], xlabel="time_SIM [s]", ylabel="BOFq_SSB", label=["q0","q1","q2","q3"], title="Spacecraft Attitude")

def DYN_STR_plot(timeline, DYN_obj, SEN_obj, NAV_obj):
    STARSdir_SSB = timeline["DYN"]["DYN_STR"]["STARSdir_SSB"][0] # [time, star, direction]

    # Reshape stars into a big list of direction vectors
    STARSdir_SSB_reshaped = STARSdir_SSB.reshape(-1, 3) # [time * star, direction]
    x = STARSdir_SSB_reshaped[:,0]
    y = STARSdir_SSB_reshaped[:,1]
    z = STARSdir_SSB_reshaped[:,2]

    # 3D sky sphere
    PPC.plot(x, y, z, style='.', label="Stars", xlabel="X SSB", ylabel="Y SSB", zlabel="Z SSB", title="Star Field Normalized", aspect="equal", color=colors["black"])

    # 2D sky sphere
    STARSproj = PPC.aitoff_projection(STARSdir_SSB_reshaped)
    boundary  = PPC.aitoff_boundary()
    fig, ax = PPC.plot(STARSproj[0], STARSproj[1], style='.', label="Stars", xlabel="Right Ascension [deg]", ylabel="Declination [deg]", title="Star Field - Aitoff Projection", aspect="equal", color=colors["black"])
    PPC.plot(boundary[0], boundary[1], fig=fig, ax=ax)

def DYN_PSR_plot(timeline, DYN_obj, SEN_obj, NAV_obj):
    time_SIM       = timeline["DYN"]["DYN_TIME"]["time_SIM"]
    STARSdir_SSB   = timeline["DYN"]["DYN_STR"]["STARSdir_SSB"][0] # [time, star, direction]
    PULSARSdir_SSB = timeline["DYN"]["DYN_PSR"]["PULSARSdir_SSB"][0] # [time, pulsar, direction]
    roemer_delay   = timeline["DYN"]["DYN_PSR"]["roemer_delay"]
    shapiro_delay  = timeline["DYN"]["DYN_PSR"]["shapiro_delay"]
    SCdt_SSB       = timeline["DYN"]["DYN_PSR"]["SCdt_SSB"]
    PULSARname     = DYN_obj.DYN_PSR.par["name"]
    n_pulsars      = DYN_obj.DYN_PSR.par["n_pulsars"]

    # Reshape stars into a big list of direction vectors
    STARSdir_SSB_reshaped = STARSdir_SSB.reshape(-1, 3) # [time * star, direction]
    x_star = STARSdir_SSB_reshaped[:,0]
    y_star = STARSdir_SSB_reshaped[:,1]
    z_star = STARSdir_SSB_reshaped[:,2]

    # Reshape pulsars into a big list of direction vectors
    PULSARSdir_SSB_reshaped = PULSARSdir_SSB.reshape(-1, 3) # [time * star, direction]
    x_pulsar = PULSARSdir_SSB_reshaped[:,0]
    y_pulsar = PULSARSdir_SSB_reshaped[:,1]
    z_pulsar = PULSARSdir_SSB_reshaped[:,2]

    # 3D sky sphere
    fig, ax = PPC.plot(x_star, y_star, z_star, style='.', label="Stars",   xlabel="X SSB", ylabel="Y SSB", zlabel="Z SSB", title="Star Field Normalized", aspect="equal", color=colors["black"])
    PPC.plot(x_pulsar, y_pulsar, z_pulsar, style='x', label="Pulsars", color=colors["purple"], fig=fig, ax=ax)

    # 2D sky sphere
    STARSproj   = PPC.aitoff_projection(STARSdir_SSB_reshaped)
    PULSARSproj = PPC.aitoff_projection(PULSARSdir_SSB_reshaped)
    boundary    = PPC.aitoff_boundary()
    fig, ax = PPC.plot(STARSproj[0], STARSproj[1], style='.', label="Stars", xlabel="Right Ascension [deg]", ylabel="Declination [deg]", title="Star Field - Aitoff Projection", aspect="equal", color=colors["black"])
    PPC.plot(PULSARSproj[0], PULSARSproj[1], style='x', label="Pulsars", color=colors["purple"], fig=fig, ax=ax)
    PPC.plot(boundary[0], boundary[1], fig=fig, ax=ax)

    # Delays for TOAs between SC and SSB
    fig, ax1 = PPC.plot([], [], ylabel="SCdt_SSB [s]", title="True TOA delay on SC", subplot=(3,1,1))
    fig, ax2 = PPC.plot([], [], ylabel="roemer_delay [s]", fig=fig, subplot=(3,1,2))
    fig, ax3 = PPC.plot([], [], xlabel="time_SIM [s]", ylabel="shapiro_delay [s]", fig=fig, subplot=(3,1,3))
    for i in range(n_pulsars):
        PPC.plot(time_SIM, SCdt_SSB[:,i], label=PULSARname[i], xlabel="time_SIM [s]", ylabel="SCdt_SSB [s]", title="True TOA delay on SC", fig=fig, ax=ax1)
        PPC.plot(time_SIM, roemer_delay[:,i], label=PULSARname[i], ylabel="roemer_delay [s]", fig=fig, ax=ax2)
        PPC.plot(time_SIM, shapiro_delay[:,i], label=PULSARname[i], ylabel="shapiro_delay [s]", fig=fig, ax=ax3)

def SEN_STR_plot(timeline, DYN_obj, SEN_obj, NAV_obj):
    time_SIM     = timeline["DYN"]["DYN_TIME"]["time_SIM"]
    time_TDB     = timeline["DYN"]["DYN_TIME"]["time_TDB"]
    BOFq_SSB     = timeline["DYN"]["DYN_ATT"]["BOFq_SSB"]
    STARSdir_SSB = timeline["DYN"]["DYN_STR"]["STARSdir_SSB"][0] # [time, star, direction]

    time_STR      = timeline["SEN"]["SEN_STR"]["time_STR"]
    STRoutflg     = timeline["SEN"]["SEN_STR"]["STRoutflg"]
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
    fig, ax = PPC.plot(time_SIM, STRoutflg, xlabel="time_SIM [s]", ylabel="STRoutflg", label="STRoutflag", title="STR output flag")

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

    fig, ax = PPC.plot([], [], xlabel="X STR", ylabel="Y STR", style='.', label="Stars", aspect="equal", color=colors["green"])
    PPC.plot(boundary[0], boundary[1], label="FOV", style='--', fig=fig, ax=ax, color=colors["black"])
    if not PPC.is_nan(SUNproj):    PPC.plot(SUNproj[0],    SUNproj[1],    label="Sun",    style='.',  fig=fig, ax=ax, color=colors["orange"])
    if not PPC.is_nan(MOONproj):   PPC.plot(MOONproj[0],   MOONproj[1],   label="Moon",   style='.',  fig=fig, ax=ax, color=colors["grey"])
    if not PPC.is_nan(EARTHproj):  PPC.plot(EARTHproj[0],  EARTHproj[1],  label="Earth",  style='.',  fig=fig, ax=ax, color=colors["blue"])
    if not PPC.is_nan(DEIMOSproj): PPC.plot(DEIMOSproj[0], DEIMOSproj[1], label="Deimos", style='.',  fig=fig, ax=ax, color=colors["lightgrey"])
    if not PPC.is_nan(PHOBOSproj): PPC.plot(PHOBOSproj[0], PHOBOSproj[1], label="Phobos", style='.',  fig=fig, ax=ax, color=colors["darkgrey"])
    if not PPC.is_nan(MARSproj):   PPC.plot(MARSproj[0],   MARSproj[1],   label="Mars",   style='.',  fig=fig, ax=ax, color=colors["red"])
    if not PPC.is_nan(STARSproj):  PPC.plot(STARSproj[0],  STARSproj[1],                  style='.',  fig=fig, ax=ax, color=colors["green"])

    # 3D sky sphere
    fig, ax = PPC.plot(STARSdir_SSB[:,0], STARSdir_SSB[:,1], STARSdir_SSB[:,2], style='.', label="Stars", xlabel="X SC", ylabel="Y SC", zlabel="Z SC", title="Star Field Normalized in SC frame", aspect="equal", color=colors["black"])
    if not PPC.is_nan(STARSdir_SC_mes):   PPC.plot(STARSdir_SC_mes[:,0],  STARSdir_SC_mes[:,1],  STARSdir_SC_mes[:,2],  style='.', label="Visible Stars", fig=fig, ax=ax, color=colors["green"])
    if not PPC.is_nan(SUNdir_SC_mes):     PPC.plot(SUNdir_SC_mes[:,0],    SUNdir_SC_mes[:,1],    SUNdir_SC_mes[:,1],    style='.', label="Sun",           fig=fig, ax=ax, color=colors["orange"])
    if not PPC.is_nan(MOONdir_SC_mes):    PPC.plot(MOONdir_SC_mes[:,0],   MOONdir_SC_mes[:,1],   MOONdir_SC_mes[:,2],   style='.', label="MOON",          fig=fig, ax=ax, color=colors["grey"])
    if not PPC.is_nan(EARTHdir_SC_mes):   PPC.plot(EARTHdir_SC_mes[:,0],  EARTHdir_SC_mes[:,1],  EARTHdir_SC_mes[:,2],  style='.', label="Earth",         fig=fig, ax=ax, color=colors["blue"])
    if not PPC.is_nan(DEIMOSdir_SC_mes):  PPC.plot(DEIMOSdir_SC_mes[:,0], DEIMOSdir_SC_mes[:,1], DEIMOSdir_SC_mes[:,2], style='.', label="Deimos",        fig=fig, ax=ax, color=colors["lightgrey"])
    if not PPC.is_nan(PHOBOSdir_SC_mes):  PPC.plot(PHOBOSdir_SC_mes[:,0], PHOBOSdir_SC_mes[:,1], PHOBOSdir_SC_mes[:,2], style='.', label="Phobos",        fig=fig, ax=ax, color=colors["darkgrey"])
    if not PPC.is_nan(MARSdir_SC_mes):    PPC.plot(MARSdir_SC_mes[:,0],   MARSdir_SC_mes[:,1],   MARSdir_SC_mes[:,2],   style='.', label="Mars",          fig=fig, ax=ax, color=colors["red"])
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

def SEN_PSR_plot(timeline, DYN_obj, SEN_obj, NAV_obj):
    time_SIM = timeline["DYN"]["DYN_TIME"]["time_SIM"]
    SCdt_SSB = timeline["DYN"]["DYN_PSR"]["SCdt_SSB"]

    # TODO: Plot the other SEN_PSR outputs
    time_PSR     = timeline["SEN"]["SEN_PSR"]["time_PSR"]
    PSRoutflg    = timeline["SEN"]["SEN_PSR"]["PSRoutflg"]
    SCdt_SSB_mes = timeline["SEN"]["SEN_PSR"]["SCdt_SSB_mes"]

    PULSARname = DYN_obj.DYN_PSR.par["name"]
    n_pulsars  = SEN_obj.SEN_PSR.par["n_pulsars"]

    light_speed_cst = CONSTANTS_par["light_speed_cst"]

    # Pulsar signals on the SC
    PULSARname = DYN_obj.DYN_PSR.par["name"]

    # Measured vs True delay for TOAs between SC and SSB
    fig, ax1 = PPC.plot([], [], ylabel="SCdt_SSB_mes [s]", title="Measured vs True TOA delay on SC", subplot=(3,1,1))
    fig, ax2 = PPC.plot([], [], ylabel="Time noise [s]", fig=fig, subplot=(3,1,2))
    fig, ax3 = PPC.plot([], [], xlabel="time_SIM [s]", ylabel="Range noise [km]", fig=fig, subplot=(3,1,3))
    time_noise = SCdt_SSB - SCdt_SSB_mes
    range_noise = light_speed_cst*time_noise
    for i in range(n_pulsars):
        PPC.plot(time_SIM, SCdt_SSB_mes[:,i], label=f"{PULSARname[i]} Meas",        fig=fig, ax=ax1)
        PPC.plot(time_SIM, time_noise[:,i],   label=f"{PULSARname[i]} time noise",  fig=fig, ax=ax2)
        PPC.plot(time_SIM, range_noise[:,i],  label=f"{PULSARname[i]} range noise", fig=fig, ax=ax3)

def NAV_EPH_plot(timeline, DYN_obj, SEN_obj, NAV_obj):
    time_SIM = timeline["DYN"]["DYN_TIME"]["time_SIM"]

    bodies = [
        dict(name="SUN",    idx=0, color=colors["orange"]),
        dict(name="MOON",   idx=2, color=colors["grey"]),
        dict(name="EARTH",  idx=1, color=colors["blue"]),
        dict(name="DEIMOS", idx=4, color=colors["darkgrey"]),
        dict(name="PHOBOS", idx=5, color=colors["lightgrey"]),
        dict(name="MARS",   idx=3, color=colors["red"]),
    ]

    for body in bodies:
        name = body["name"]
        # Plot reference position from DYN
        fig, ax1 = PPC.plot(time_SIM, timeline["DYN"]["DYN_EPH"][f"{name}pos_SSB"], ylabel=f"{name}pos_SSB [km]", label=["DYN x","DYN y","DYN z"], title=f"{name}pos_SSB", style='.', subplot=(2,1,1))
        # Plot position from NAV_EPH
        PPC.plot(time_SIM, timeline["NAV"]["NAV_EPH"][f"{name}pos_SSB"], xlabel="time_SIM [s]", ylabel=f"{name}pos_SSB [km]", label=["NAV x","NAV y","NAV z"], title=f"{name}pos_SSB", fig=fig, ax=ax1)
        # Plot reference velocity from DYN
        fig, ax2 = PPC.plot(time_SIM, timeline["DYN"]["DYN_EPH"][f"{name}vel_SSB"], ylabel=f"{name}vel_SSB [km]", label=["DYN Vx","DYN vy","DYN vz"], title=f"{name}vel_SSB", style='.', fig=fig, subplot=(2,1,2))
        # Plot velocity from NAV_EPH
        PPC.plot(time_SIM, timeline["NAV"]["NAV_EPH"][f"{name}vel_SSB"], xlabel="time_SIM [s]", ylabel=f"{name}vel_SSB [km/s]", label=["NAV x","NAV y","NAV z"], title=f"{name}vel_SSB", fig=fig, ax=ax2)

def NAV_STR_plot(timeline, DYN_obj, SEN_obj, NAV_obj):
    DYN_STARSdir_SSB = timeline["DYN"]["DYN_STR"]["STARSdir_SSB"][0] # [time, star, direction]
    NAV_STARSdir_SSB = timeline["NAV"]["NAV_STR"]["STARSdir_SSB"][0] # [time, star, direction]

    # Reshape stars into a big list of direction vectors
    DYN_STARSdir_SSB_reshaped = DYN_STARSdir_SSB.reshape(-1, 3) # [time * star, direction]
    NAV_STARSdir_SSB_reshaped = NAV_STARSdir_SSB.reshape(-1, 3) # [time * star, direction]
    DYN_x = DYN_STARSdir_SSB_reshaped[:,0]
    DYN_y = DYN_STARSdir_SSB_reshaped[:,1]
    DYN_z = DYN_STARSdir_SSB_reshaped[:,2]
    NAV_x = NAV_STARSdir_SSB_reshaped[:,0]
    NAV_y = NAV_STARSdir_SSB_reshaped[:,1]
    NAV_z = NAV_STARSdir_SSB_reshaped[:,2]

    # 3D sky sphere
    fig, ax = PPC.plot(DYN_x, DYN_y, DYN_z, style='.', label="Stars DYN", xlabel="X SSB", ylabel="Y SSB", zlabel="Z SSB", title="Star Field Normalized", aspect="equal", color=colors["black"])
    PPC.plot(NAV_x, NAV_y, NAV_z, style='.', label="Stars NAV", color=colors["blue"], fig=fig, ax=ax)

    # 2D sky sphere
    DYN_STARSproj = PPC.aitoff_projection(DYN_STARSdir_SSB_reshaped)
    NAV_STARSproj = PPC.aitoff_projection(NAV_STARSdir_SSB_reshaped)
    boundary  = PPC.aitoff_boundary()
    fig, ax = PPC.plot(DYN_STARSproj[0], DYN_STARSproj[1], style='.', label="DYN Stars", xlabel="Right Ascension [deg]", ylabel="Declination [deg]", title="Star Field - Aitoff Projection", aspect="equal", color=colors["black"])
    fig, ax = PPC.plot(NAV_STARSproj[0], NAV_STARSproj[1], style='.', label="NAV Stars", color=colors["blue"], fig=fig, ax=ax)
    PPC.plot(boundary[0], boundary[1], fig=fig, ax=ax)

def NAV_CEL_plot(timeline, DYN_obj, SEN_obj, NAV_obj):
    time_SIM        = timeline["DYN"]["DYN_TIME"]["time_SIM"]
    SCpos_SSB       = timeline["DYN"]["DYN_TRA"]["SCpos_SSB"]
    SCvel_SSB       = timeline["DYN"]["DYN_TRA"]["SCvel_SSB"]
    STRoutflg       = timeline["SEN"]["SEN_STR"]["STRoutflg"]
    STARSdir_SC_mes = timeline["SEN"]["SEN_STR"]["STARSdir_SC_mes"]

    NAV_CELoutflg               = timeline["NAV"]["NAV_CEL"]["NAV_CELoutflg"]
    BODYsel_STARdir_SC_mes_list = timeline["NAV"]["NAV_CEL"]["BODYsel_STARdir_SC_mes_list"]

    z = timeline["NAV"]["NAV_CEL"]["z"]
    R = timeline["NAV"]["NAV_CEL"]["R"]

    bodies = [
        dict(name="SUN",    idx=0, color=colors["orange"]),
        dict(name="MOON",   idx=2, color=colors["grey"]),
        dict(name="EARTH",  idx=1, color=colors["blue"]),
        dict(name="DEIMOS", idx=4, color=colors["darkgrey"]),
        dict(name="PHOBOS", idx=5, color=colors["lightgrey"]),
        dict(name="MARS",   idx=3, color=colors["red"]),
    ]

    angles_deg = {}

    for body in bodies:
        idx = body["idx"]
        angles_deg[body["name"]] = (np.cos(z[:, idx])*CONSTANTS_par["rad2deg_cst"])

    active_bodies = [
        body for body in bodies if not PPC.is_nan(angles_deg[body["name"]])
    ]

    # Reshape vectors for plot
    STARSdir_SC_mes_reshaped        = STARSdir_SC_mes.reshape(-1, 3) # [time * star, direction]

    # Plot status
    fig, ax = PPC.plot(time_SIM, NAV_CELoutflg, xlabel="time_SIM [s]", ylabel="flag", label="NAV_CELoutflag", title="NAV_CEL output flag")
    PPC.plot(time_SIM, STRoutflg, label="STRoutflag", style='--', fig=fig, ax=ax)

    # Plot angles (arccosine of measurement model)
    fig = None
    ax  = None
    for body in active_bodies:
        name = body["name"]
        angle = angles_deg[name]
        fig, ax = PPC.plot(time_SIM, angle, xlabel="time_SIM [s]", ylabel=f"{name}angles_mes [deg]", label=f"Angle {name}", title=f"Line-of-Sight Angle from {name} to Star", color=body["color"], fig=fig, ax=ax)

    # Plot covariance
    R_diag  = np.diagonal(R, axis1=1, axis2=2)
    sigma_z = np.sqrt(R_diag)
    fig = None
    for body in active_bodies:
        idx  = body["idx"]
        name = body["name"]
        fig, ax = PPC.plot(time_SIM, sigma_z[:, idx], label=f"σ_z {name}", xlabel="time_SIM [s]", ylabel="σ", title="Measurement Standard Deviation", color=body["color"], fig=fig)

    # Plot predicted measurement, assuming the true spacecraft position as the state + innovation
    x_true   = np.hstack((SCpos_SSB, SCvel_SSB))
    n_iter   = len(x_true)
    n_bodies = len(bodies)
    h_hist   = np.zeros((n_iter, n_bodies))
    H_hist   = np.zeros((n_iter, n_bodies, 6))
    for k in range(n_iter):
        # First update NAV_CEL state, otherwise h(x) will be computed for the last (already computed) state
        NAV_obj.NAV_CEL.state["BODYpos_SSB_list"] = timeline["NAV"]["NAV_CEL"]["BODYpos_SSB_list"][k]
        NAV_obj.NAV_CEL.state["BODYsel_STARdir_SC_ref_list"] = timeline["NAV"]["NAV_CEL"]["BODYsel_STARdir_SC_ref_list"][k]
        h_hist[k] = NAV_obj.NAV_CEL.h(x_true[k])
        H_hist[k] = NAV_obj.NAV_CEL.H(x_true[k])
    fig = None
    ax1 = None
    ax2 = None
    innov = z - h_hist
    for body in active_bodies:
        name = body["name"]
        idx  = body["idx"]
        fig, ax1 = PPC.plot(time_SIM, h_hist[:, idx], label=f"h(x) {name}", ylabel="h(x)", title="Predicted Measurement h(x)", color=body["color"], fig=fig, ax=ax1, subplot=(2,1,1))
        PPC.plot(time_SIM, z[:, idx], label=f"z {name}", style='--', color=body["color"], fig=fig, ax=ax1)
        fig, ax2 = PPC.plot(time_SIM, innov[:, idx], label=f"{name}", xlabel="time_SIM [s]", ylabel="z - h(x)", title="Innovation z - h(x)", color=body["color"], fig=fig, ax=ax2, subplot=(2,1,2))

    # Plot Selected stars for each body
    fig, ax = PPC.plot(STARSdir_SC_mes_reshaped[:,0], STARSdir_SC_mes_reshaped[:,1], STARSdir_SC_mes_reshaped[:,2], style='.', label="Visible Stars", xlabel="X SC", ylabel="Y SC", zlabel="Z SC", title="Star Field Normalized in SC frame", aspect="equal", color=colors["green"])
    PPC.plot([0], [0], [0],  style='+', label="SC",  fig=fig, ax=ax, color=colors["magenta"])

    for body in active_bodies:
        name = body["name"]
        idx  = body["idx"]
        BODYsel_STARdir_SC_mes = BODYsel_STARdir_SC_mes_list[:,idx,:] # [time, body, direction]
        BODYdir_SC_mes         = timeline["SEN"]["SEN_STR"][f"{name}dir_SC_mes"]

        PPC.plot(BODYsel_STARdir_SC_mes[:,0], BODYsel_STARdir_SC_mes[:,1], BODYsel_STARdir_SC_mes[:,2], style='.', label=f"{name} selected star", fig=fig, ax=ax, color=body["color"])
        PPC.plot(BODYdir_SC_mes[:,0],         BODYdir_SC_mes[:,1],         BODYdir_SC_mes[:,2],         style='x', label=f"{name}",               fig=fig, ax=ax, color=body["color"])

    # Plot selected stars direction over time
    fig, ax1 = PPC.plot(time_SIM, STARSdir_SC_mes[:, :, 0], color=colors["green"], subplot=(3,1,1), ylabel="x", title=f"Direction of Visible Stars and Reference Selected Star")
    fig, ax2 = PPC.plot(time_SIM, STARSdir_SC_mes[:, :, 1], color=colors["green"], subplot=(3,1,2), ylabel="y", fig=fig)
    fig, ax3 = PPC.plot(time_SIM, STARSdir_SC_mes[:, :, 2], color=colors["green"], subplot=(3,1,3), xlabel="time_SIM [s]", ylabel="z", fig=fig)
    for body in active_bodies:
        name = body["name"]
        idx  = body["idx"]
        BODYsel_STARdir_SC_mes = BODYsel_STARdir_SC_mes_list[:,idx,:] # [time, body, direction]
        PPC.plot(time_SIM, BODYsel_STARdir_SC_mes[:, 0], label=f"{name}", color=body["color"], fig=fig, ax=ax1)
        PPC.plot(time_SIM, BODYsel_STARdir_SC_mes[:, 1], label=f"{name}", color=body["color"], fig=fig, ax=ax2)
        PPC.plot(time_SIM, BODYsel_STARdir_SC_mes[:, 2], label=f"{name}", color=body["color"], fig=fig, ax=ax3)

PPC_plots = {
    "DYN_TIME"  : DYN_TIME_plot,
    "DYN_EPH"   : DYN_EPH_plot,
    "DYN_GRV"   : DYN_GRV_plot,
    "DYN_TRA"   : DYN_TRA_plot,
    "DYN_ATT"   : DYN_ATT_plot,
    "DYN_STR"   : DYN_STR_plot,
    "DYN_PSR"   : DYN_PSR_plot,
    "SEN_STR"   : SEN_STR_plot,
    "SEN_PSR"   : SEN_PSR_plot,
    "NAV_EPH"   : NAV_EPH_plot,
    "NAV_STR"   : NAV_STR_plot,
    "NAV_CEL"   : NAV_CEL_plot,
}
