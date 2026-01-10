# List all plots to be generated, organized by group name, to be selected on SIM_par

from PPC import PPC

# Define color cycle
colors = {
    "blue"      : "#0066FF",
    "red"       : "#CC0000",
    "green"     : "#33CC00",
    "magenta"   : "#FF00FF",
    "orange"    : "#FE9920",
    "grey"      : "#505050",
    "lightgrey" : "#707070",
    "darkgrey"  : "#303030",
    "black"     : "#111111",
}
PPC.setup_plot(colors.values())

def DYN_TIME_plot(timeline):
    PPC.plot(timeline["DYN"]["DYN_TIME"]["time_SIM"], timeline["DYN"]["DYN_TIME"]["time_TDB"], xlabel="time_SIM [s]", ylabel="time_TDB [s]", title="Time")

def DYN_SUN_plot(timeline):
    # Sun state
    fig, ax = PPC.plot(timeline["DYN"]["DYN_TIME"]["time_SIM"], timeline["DYN"]["DYN_SUN"]["SUNpos_SSB"], ylabel="SUNpos_SSB [km]", label=["x","y","z"], title="SUNpos_SSB", subplot=(2,1,1))
    PPC.plot(timeline["DYN"]["DYN_TIME"]["time_SIM"], timeline["DYN"]["DYN_SUN"]["SUNvel_SSB"], xlabel="time_SIM [s]", ylabel="SUNvel_SSB [km]", label=["x","y","z"], title="SUNvel_SSB", fig=fig, subplot=(2,1,2))

def DYN_EARTH_plot(timeline):
    # Earth state
    fig, ax = PPC.plot(timeline["DYN"]["DYN_TIME"]["time_SIM"], timeline["DYN"]["DYN_EARTH"]["EARTHpos_SSB"], ylabel="EARTHpos_SSB [km]", label=["x","y","z"], title="EARTHpos_SSB", subplot=(2,1,1))
    PPC.plot(timeline["DYN"]["DYN_TIME"]["time_SIM"], timeline["DYN"]["DYN_EARTH"]["EARTHvel_SSB"], xlabel="time_SIM [s]", ylabel="EARTHvel_SSB [km]", label=["x","y","z"], title="EARTHvel_SSB", fig=fig, subplot=(2,1,2))

    # Moon state
    fig, ax = PPC.plot(timeline["DYN"]["DYN_TIME"]["time_SIM"], timeline["DYN"]["DYN_EARTH"]["MOONpos_SSB"], ylabel="MOONpos_SSB [km]", label=["x","y","z"], title="MOONpos_SSB", subplot=(2,1,1))
    PPC.plot(timeline["DYN"]["DYN_TIME"]["time_SIM"], timeline["DYN"]["DYN_EARTH"]["MOONvel_SSB"], xlabel="time_SIM [s]", ylabel="MOONvel_SSB [km]", label=["x","y","z"], title="MOONvel_SSB", fig=fig, subplot=(2,1,2))

def DYN_MARS_plot(timeline):
    # Mars state
    fig, ax = PPC.plot(timeline["DYN"]["DYN_TIME"]["time_SIM"], timeline["DYN"]["DYN_MARS"]["MARSpos_SSB"], ylabel="MARSpos_SSB [km]", label=["x","y","z"], title="MARSpos_SSB", subplot=(2,1,1))
    PPC.plot(timeline["DYN"]["DYN_TIME"]["time_SIM"], timeline["DYN"]["DYN_MARS"]["MARSvel_SSB"], xlabel="time_SIM [s]", ylabel="MARSvel_SSB [km]", label=["x","y","z"], title="MARSvel_SSB", fig=fig, subplot=(2,1,2))

    # Deimos state
    fig, ax = PPC.plot(timeline["DYN"]["DYN_TIME"]["time_SIM"], timeline["DYN"]["DYN_MARS"]["DEIMOSpos_SSB"], ylabel="DEIMOSpos_SSB [km]", label=["x","y","z"], title="DEIMOSpos_SSB", subplot=(2,1,1))
    PPC.plot(timeline["DYN"]["DYN_TIME"]["time_SIM"], timeline["DYN"]["DYN_MARS"]["DEIMOSvel_SSB"], xlabel="time_SIM [s]", ylabel="DEIMOSvel_SSB [km]", label=["x","y","z"], title="DEIMOSvel_SSB", fig=fig, subplot=(2,1,2))

    # Phobos state
    fig, ax = PPC.plot(timeline["DYN"]["DYN_TIME"]["time_SIM"], timeline["DYN"]["DYN_MARS"]["PHOBOSpos_SSB"], ylabel="PHOBOSpos_SSB [km]", label=["x","y","z"], title="PHOBOSpos_SSB", subplot=(2,1,1))
    PPC.plot(timeline["DYN"]["DYN_TIME"]["time_SIM"], timeline["DYN"]["DYN_MARS"]["PHOBOSvel_SSB"], xlabel="time_SIM [s]", ylabel="PHOBOSvel_SSB [km]", label=["x","y","z"], title="PHOBOSvel_SSB", fig=fig, subplot=(2,1,2))

def DYN_GRV_plot(timeline):
    # Gravity acceleration on SSB frame
    PPC.plot(timeline["DYN"]["DYN_TIME"]["time_SIM"], timeline["DYN"]["DYN_GRV"]["grvacc_SSB"], ylabel="grvacc_SSB [km/s^2]", label=["x","y","z"], title="grvacc_SSB")

def DYN_TRA_plot(timeline):
    # 3D trajectories
    fig, ax = PPC.plot([0], [0], [0], style='+', label="SSB", xlabel="X SSB [km]", ylabel="Y SSB [km]", zlabel="Z SSB [km]", title="Simulation Trajectories", aspect="equal",                        color=colors["black"])
    PPC.plot(timeline["DYN"]["DYN_SUN"]["SUNpos_SSB"][:,0],     timeline["DYN"]["DYN_SUN"]["SUNpos_SSB"][:,1],     timeline["DYN"]["DYN_SUN"]["SUNpos_SSB"][:,2],     label="Sun",   fig=fig, ax=ax, color=colors["orange"])
    PPC.plot(timeline["DYN"]["DYN_EARTH"]["MOONpos_SSB"][:,0],  timeline["DYN"]["DYN_EARTH"]["MOONpos_SSB"][:,1],  timeline["DYN"]["DYN_EARTH"]["MOONpos_SSB"][:,2],  label="Moon",  fig=fig, ax=ax, color=colors["grey"])
    PPC.plot(timeline["DYN"]["DYN_EARTH"]["EARTHpos_SSB"][:,0], timeline["DYN"]["DYN_EARTH"]["EARTHpos_SSB"][:,1], timeline["DYN"]["DYN_EARTH"]["EARTHpos_SSB"][:,2], label="Earth", fig=fig, ax=ax, color=colors["blue"])
    PPC.plot(timeline["DYN"]["DYN_MARS"]["DEIMOSpos_SSB"][:,0], timeline["DYN"]["DYN_MARS"]["DEIMOSpos_SSB"][:,1], timeline["DYN"]["DYN_MARS"]["DEIMOSpos_SSB"][:,2], label="Deimos",fig=fig, ax=ax, color=colors["lightgrey"])
    PPC.plot(timeline["DYN"]["DYN_MARS"]["PHOBOSpos_SSB"][:,0], timeline["DYN"]["DYN_MARS"]["PHOBOSpos_SSB"][:,1], timeline["DYN"]["DYN_MARS"]["PHOBOSpos_SSB"][:,2], label="Phobos",fig=fig, ax=ax, color=colors["darkgrey"])
    PPC.plot(timeline["DYN"]["DYN_MARS"]["MARSpos_SSB"][:,0],   timeline["DYN"]["DYN_MARS"]["MARSpos_SSB"][:,1],   timeline["DYN"]["DYN_MARS"]["MARSpos_SSB"][:,2],   label="Mars",  fig=fig, ax=ax, color=colors["red"])
    PPC.plot(timeline["DYN"]["DYN_TRA"]["SCpos_SSB"][:,0],      timeline["DYN"]["DYN_TRA"]["SCpos_SSB"][:,1],      timeline["DYN"]["DYN_TRA"]["SCpos_SSB"][:,2],      label="SC",    fig=fig, ax=ax, color=colors["magenta"])

    # 2D SC and Moon trajectories around Earth
    MOONpos_ECI = timeline["DYN"]["DYN_EARTH"]["EARTHpos_SSB"] - timeline["DYN"]["DYN_EARTH"]["MOONpos_SSB"]
    fig, ax = PPC.plot([0], [0], style='o', label="Earth", xlabel="X ECI [km]", ylabel="Y ECI [km]", title="Trajectories around Earth", aspect="equal", color=colors["blue"])
    PPC.plot(timeline["DYN"]["DYN_TRA"]["SCpos_ECI"][:,0], timeline["DYN"]["DYN_TRA"]["SCpos_ECI"][:,1], label="SC",   fig=fig, ax=ax, color=colors["magenta"])
    PPC.plot(MOONpos_ECI[:,0],                             MOONpos_ECI[:,1],                             label="Moon", fig=fig, ax=ax, color=colors["grey"])

    # Spacecraft state
    fig, ax = PPC.plot(timeline["DYN"]["DYN_TIME"]["time_SIM"], timeline["DYN"]["DYN_TRA"]["SCpos_SSB"], ylabel="SCpos_SSB [km]", label=["x","y","z"], title="SCpos_SSB", subplot=(2,1,1))
    PPC.plot(timeline["DYN"]["DYN_TIME"]["time_SIM"], timeline["DYN"]["DYN_TRA"]["SCvel_SSB"], ylabel="SCvel_SSB [km]", label=["x","y","z"], title="SCvel_SSB", fig=fig, subplot=(2,1,2))

def DYN_ATT_plot(timeline):
    # Spacecraft attitude
    PPC.plot(timeline["DYN"]["DYN_TIME"]["time_SIM"], timeline["DYN"]["DYN_ATT"]["BOFq_SSB"], xlabel="time_SIM [s]", ylabel="BOFq_SSB", label=["q0","q1","q2","q3"], title="Spacecraft Attitude")

def DYN_STR_plot(timeline):
    # 3D sky sphere
    x, y, z = timeline["DYN"]["DYN_STR"]["STARSdir_SSB"][0].T
    PPC.plot(x, y, z, style='.', label="Stars", xlabel="X SSB", ylabel="Y SSB", zlabel="Z SSB", title="Star Field Normalized", aspect='equal', color=colors["black"])

def SEN_STR_plot(timeline):
    time_SIM = timeline["DYN"]["DYN_TIME"]["time_SIM"]
    time_STR = timeline["SEN"]["SEN_STR"]["time_STR"]

    STRoutflg = timeline["SEN"]["SEN_STR"]["STRoutflg"]

    SUNdir_STR_mes    = timeline["SEN"]["SEN_STR"]["SUNdir_STR_mes"].T
    EARTHdir_STR_mes  = timeline["SEN"]["SEN_STR"]["EARTHdir_STR_mes"].T
    MARSdir_STR_mes   = timeline["SEN"]["SEN_STR"]["MARSdir_STR_mes"].T
    STARSdir_STR_mes  = timeline["SEN"]["SEN_STR"]["STARSdir_STR_mes"].T # [coord, star, time]

    SUNdir_SSB_mes    = timeline["SEN"]["SEN_STR"]["SUNdir_SSB_mes"].T
    EARTHdir_SSB_mes  = timeline["SEN"]["SEN_STR"]["EARTHdir_SSB_mes"].T
    MARSdir_SSB_mes   = timeline["SEN"]["SEN_STR"]["MARSdir_SSB_mes"].T
    STARSdir_SSB_mes  = timeline["SEN"]["SEN_STR"]["STARSdir_SSB_mes"].T

    STARSdir_SSB = timeline["DYN"]["DYN_STR"]["STARSdir_SSB"].T

    # Compare STR and SIM times
    fig, ax = PPC.plot(time_SIM, time_SIM, xlabel="time_SIM [s]", ylabel="time [s]", label="time_SIM", title="STR Time")
    PPC.plot(time_SIM, time_STR, ylabel="time [s]", label="time_STR", fig=fig, ax=ax)

    # Plot STR Gnomonic lens projection
    SUNproj   = PPC.gnomonic_projection(SUNdir_STR_mes)
    EARTHproj = PPC.gnomonic_projection(EARTHdir_STR_mes)
    MARSproj  = PPC.gnomonic_projection(MARSdir_STR_mes)
    STARSproj = PPC.gnomonic_projection(STARSdir_STR_mes)

    fig, ax = PPC.plot([], [], xlabel="X STR", ylabel="Y STR", label="Stars", aspect="equal", color=colors["green"])
    PPC.plot(SUNproj[0],   SUNproj[1],    label="Sun",   style='.', fig=fig, ax=ax, color=colors["orange"])
    PPC.plot(EARTHproj[0], EARTHproj[1],  label="Earth", style='.', fig=fig, ax=ax, color=colors["blue"])
    PPC.plot(MARSproj[0],  MARSproj[1],   label="Mars",  style='.', fig=fig, ax=ax, color=colors["red"])
    PPC.plot(STARSproj[0], STARSproj[1],                 style='.', fig=fig, ax=ax, color=colors["green"])

    # 3D sky sphere on step 0
    fig, ax = PPC.plot(STARSdir_SSB[0,:,1], STARSdir_SSB[1,:,1], STARSdir_SSB[2,:,1], style='.', label="Stars", xlabel="X SSB", ylabel="Y SSB", zlabel="Z SSB", title="Star Field Normalized in SSB frame", aspect='equal', color=colors["black"])
    PPC.plot(STARSdir_SSB_mes[0,:,1], STARSdir_SSB_mes[1,:,1], STARSdir_SSB_mes[2,:,1],   style='.', label="Visible Stars", fig=fig, ax=ax, color=colors["green"])
    PPC.plot(SUNdir_SSB_mes[0,:],     SUNdir_SSB_mes[0,:],     SUNdir_SSB_mes[0,:],   style='.-', label="Sun",   fig=fig, ax=ax, color=colors["orange"])
    PPC.plot(EARTHdir_SSB_mes[0,:],   EARTHdir_SSB_mes[0,:],   EARTHdir_SSB_mes[0,:], style='.-', label="Earth", fig=fig, ax=ax, color=colors["blue"])
    PPC.plot(MARSdir_SSB_mes[0,:],    MARSdir_SSB_mes[0,:],    MARSdir_SSB_mes[0,:],  style='.-', label="Mars",  fig=fig, ax=ax, color=colors["red"])
    PPC.plot([0], [0], [0],  style='+', label="SC",  fig=fig, ax=ax, color=colors["magenta"])

PPC_plots = {
    "DYN_TIME"  : DYN_TIME_plot,
    "DYN_SUN"   : DYN_SUN_plot,
    "DYN_EARTH" : DYN_EARTH_plot,
    "DYN_MARS"  : DYN_MARS_plot,
    "DYN_GRV"   : DYN_GRV_plot,
    "DYN_TRA"   : DYN_TRA_plot,
    "DYN_ATT"   : DYN_ATT_plot,
    "DYN_STR"   : DYN_STR_plot,
    "SEN_STR"   : SEN_STR_plot,
}
