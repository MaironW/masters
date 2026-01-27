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

def DYN_TIME_plot(timeline, DYN_obj, SEN_obj, NAV_obj):
    PPC.plot(timeline["DYN"]["DYN_TIME"]["time_SIM"], timeline["DYN"]["DYN_TIME"]["time_TDB"], xlabel="time_SIM [s]", ylabel="time_TDB [s]", title="Time")

def DYN_SUN_plot(timeline, DYN_obj, SEN_obj, NAV_obj):
    # Sun state
    fig, ax = PPC.plot(timeline["DYN"]["DYN_TIME"]["time_SIM"], timeline["DYN"]["DYN_SUN"]["SUNpos_SSB"], ylabel="SUNpos_SSB [km]", label=["x","y","z"], title="SUNpos_SSB", subplot=(2,1,1))
    PPC.plot(timeline["DYN"]["DYN_TIME"]["time_SIM"], timeline["DYN"]["DYN_SUN"]["SUNvel_SSB"], xlabel="time_SIM [s]", ylabel="SUNvel_SSB [km]", label=["x","y","z"], title="SUNvel_SSB", fig=fig, subplot=(2,1,2))

def DYN_EARTH_plot(timeline, DYN_obj, SEN_obj, NAV_obj):
    # Earth state
    fig, ax = PPC.plot(timeline["DYN"]["DYN_TIME"]["time_SIM"], timeline["DYN"]["DYN_EARTH"]["EARTHpos_SSB"], ylabel="EARTHpos_SSB [km]", label=["x","y","z"], title="EARTHpos_SSB", subplot=(2,1,1))
    PPC.plot(timeline["DYN"]["DYN_TIME"]["time_SIM"], timeline["DYN"]["DYN_EARTH"]["EARTHvel_SSB"], xlabel="time_SIM [s]", ylabel="EARTHvel_SSB [km]", label=["x","y","z"], title="EARTHvel_SSB", fig=fig, subplot=(2,1,2))

    # Moon state
    fig, ax = PPC.plot(timeline["DYN"]["DYN_TIME"]["time_SIM"], timeline["DYN"]["DYN_EARTH"]["MOONpos_SSB"], ylabel="MOONpos_SSB [km]", label=["x","y","z"], title="MOONpos_SSB", subplot=(2,1,1))
    PPC.plot(timeline["DYN"]["DYN_TIME"]["time_SIM"], timeline["DYN"]["DYN_EARTH"]["MOONvel_SSB"], xlabel="time_SIM [s]", ylabel="MOONvel_SSB [km]", label=["x","y","z"], title="MOONvel_SSB", fig=fig, subplot=(2,1,2))

def DYN_MARS_plot(timeline, DYN_obj, SEN_obj, NAV_obj):
    # Mars state
    fig, ax = PPC.plot(timeline["DYN"]["DYN_TIME"]["time_SIM"], timeline["DYN"]["DYN_MARS"]["MARSpos_SSB"], ylabel="MARSpos_SSB [km]", label=["x","y","z"], title="MARSpos_SSB", subplot=(2,1,1))
    PPC.plot(timeline["DYN"]["DYN_TIME"]["time_SIM"], timeline["DYN"]["DYN_MARS"]["MARSvel_SSB"], xlabel="time_SIM [s]", ylabel="MARSvel_SSB [km]", label=["x","y","z"], title="MARSvel_SSB", fig=fig, subplot=(2,1,2))

    # Deimos state
    fig, ax = PPC.plot(timeline["DYN"]["DYN_TIME"]["time_SIM"], timeline["DYN"]["DYN_MARS"]["DEIMOSpos_SSB"], ylabel="DEIMOSpos_SSB [km]", label=["x","y","z"], title="DEIMOSpos_SSB", subplot=(2,1,1))
    PPC.plot(timeline["DYN"]["DYN_TIME"]["time_SIM"], timeline["DYN"]["DYN_MARS"]["DEIMOSvel_SSB"], xlabel="time_SIM [s]", ylabel="DEIMOSvel_SSB [km]", label=["x","y","z"], title="DEIMOSvel_SSB", fig=fig, subplot=(2,1,2))

    # Phobos state
    fig, ax = PPC.plot(timeline["DYN"]["DYN_TIME"]["time_SIM"], timeline["DYN"]["DYN_MARS"]["PHOBOSpos_SSB"], ylabel="PHOBOSpos_SSB [km]", label=["x","y","z"], title="PHOBOSpos_SSB", subplot=(2,1,1))
    PPC.plot(timeline["DYN"]["DYN_TIME"]["time_SIM"], timeline["DYN"]["DYN_MARS"]["PHOBOSvel_SSB"], xlabel="time_SIM [s]", ylabel="PHOBOSvel_SSB [km]", label=["x","y","z"], title="PHOBOSvel_SSB", fig=fig, subplot=(2,1,2))

def DYN_GRV_plot(timeline, DYN_obj, SEN_obj, NAV_obj):
    # Gravity acceleration on SSB frame
    PPC.plot(timeline["DYN"]["DYN_TIME"]["time_SIM"], timeline["DYN"]["DYN_GRV"]["grvacc_SSB"], ylabel="grvacc_SSB [km/s^2]", label=["x","y","z"], title="grvacc_SSB")

def DYN_TRA_plot(timeline, DYN_obj, SEN_obj, NAV_obj):
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

def DYN_ATT_plot(timeline, DYN_obj, SEN_obj, NAV_obj):
    # Spacecraft attitude
    PPC.plot(timeline["DYN"]["DYN_TIME"]["time_SIM"], timeline["DYN"]["DYN_ATT"]["BOFq_SSB"], xlabel="time_SIM [s]", ylabel="BOFq_SSB", label=["q0","q1","q2","q3"], title="Spacecraft Attitude")

def DYN_STR_plot(timeline, DYN_obj, SEN_obj, NAV_obj):
    STARSdir_SSB = timeline["DYN"]["DYN_STR"]["STARSdir_SSB"][0].T

    # 3D sky sphere
    x, y, z = STARSdir_SSB
    PPC.plot(x, y, z, style='.', label="Stars", xlabel="X SSB", ylabel="Y SSB", zlabel="Z SSB", title="Star Field Normalized", aspect="equal", color=colors["black"])

    # 2D sky sphere
    STARSproj = PPC.aitoff_projection(STARSdir_SSB)
    boundary  = PPC.aitoff_boundary()
    fig, ax = PPC.plot(STARSproj[0], STARSproj[1], style='.', label="Stars", xlabel="Right Ascension [deg]", ylabel="Declination [deg]", title="Star Field - Aitoff Projection", aspect="equal", color=colors["black"])
    PPC.plot(boundary[0], boundary[1], fig=fig, ax=ax)

def SEN_STR_plot(timeline, DYN_obj, SEN_obj, NAV_obj):
    time_SIM = timeline["DYN"]["DYN_TIME"]["time_SIM"]
    time_STR = timeline["SEN"]["SEN_STR"]["time_STR"]

    STRoutflg = timeline["SEN"]["SEN_STR"]["STRoutflg"]

    SUNdir_STR_mes    = timeline["SEN"]["SEN_STR"]["SUNdir_STR_mes"].T
    EARTHdir_STR_mes  = timeline["SEN"]["SEN_STR"]["EARTHdir_STR_mes"].T
    MARSdir_STR_mes   = timeline["SEN"]["SEN_STR"]["MARSdir_STR_mes"].T
    STARSdir_STR_mes  = timeline["SEN"]["SEN_STR"]["STARSdir_STR_mes"].T # [coord, star, time]

    SUNdir_SC_mes    = timeline["SEN"]["SEN_STR"]["SUNdir_SC_mes"].T
    EARTHdir_SC_mes  = timeline["SEN"]["SEN_STR"]["EARTHdir_SC_mes"].T
    MARSdir_SC_mes   = timeline["SEN"]["SEN_STR"]["MARSdir_SC_mes"].T
    STARSdir_SC_mes  = timeline["SEN"]["SEN_STR"]["STARSdir_SC_mes"].T

    STARSdir_SSB = timeline["DYN"]["DYN_STR"]["STARSdir_SSB"].T

    field_of_view = SEN_obj.SEN_STR.par["field_of_view"]

    # Plot status
    fig, ax = PPC.plot(time_SIM, STRoutflg, xlabel="time_SIM [s]", ylabel="STRoutflg", label="STRoutflag", title="STR output flag")

    # Compare STR and SIM times
    fig, ax = PPC.plot(time_SIM, time_SIM, xlabel="time_SIM [s]", ylabel="time [s]", label="time_SIM", title="STR Time")
    PPC.plot(time_SIM, time_STR, ylabel="time [s]", label="time_STR", fig=fig, ax=ax)

    # Plot STR Gnomonic lens projection
    SUNproj   = PPC.gnomonic_projection(SUNdir_STR_mes)
    EARTHproj = PPC.gnomonic_projection(EARTHdir_STR_mes)
    MARSproj  = PPC.gnomonic_projection(MARSdir_STR_mes)
    STARSproj = PPC.gnomonic_projection(STARSdir_STR_mes)
    boundary  = PPC.gnomonic_boundary(field_of_view)

    fig, ax = PPC.plot([], [], xlabel="X STR", ylabel="Y STR", style='.', label="Stars", aspect="equal", color=colors["green"])
    PPC.plot(SUNproj[0],   SUNproj[1],   label="Sun",   style='.',  fig=fig, ax=ax, color=colors["orange"])
    PPC.plot(EARTHproj[0], EARTHproj[1], label="Earth", style='.',  fig=fig, ax=ax, color=colors["blue"])
    PPC.plot(MARSproj[0],  MARSproj[1],  label="Mars",  style='.',  fig=fig, ax=ax, color=colors["red"])
    PPC.plot(STARSproj[0], STARSproj[1],                style='.',  fig=fig, ax=ax, color=colors["green"])
    PPC.plot(boundary[0], boundary[1],   label="FOV",   style='--', fig=fig, ax=ax, color=colors["darkgrey"])

    # 3D sky sphere on step 0
    fig, ax = PPC.plot(STARSdir_SSB[0,:,1], STARSdir_SSB[1,:,1], STARSdir_SSB[2,:,1], style='.', label="Stars", xlabel="X SC", ylabel="Y SC", zlabel="Z SC", title="Star Field Normalized in SC frame", aspect="equal", color=colors["black"])
    PPC.plot(STARSdir_SC_mes[0,:,1], STARSdir_SC_mes[1,:,1], STARSdir_SC_mes[2,:,1], style='.', label="Visible Stars", fig=fig, ax=ax, color=colors["green"])
    PPC.plot(SUNdir_SC_mes[0,:],     SUNdir_SC_mes[1,:],     SUNdir_SC_mes[2,:],     style='.', label="Sun",           fig=fig, ax=ax, color=colors["orange"])
    PPC.plot(EARTHdir_SC_mes[0,:],   EARTHdir_SC_mes[1,:],   EARTHdir_SC_mes[2,:],   style='.', label="Earth",         fig=fig, ax=ax, color=colors["blue"])
    PPC.plot(MARSdir_SC_mes[0,:],    MARSdir_SC_mes[1,:],    MARSdir_SC_mes[2,:],    style='.', label="Mars",          fig=fig, ax=ax, color=colors["red"])
    PPC.plot([0], [0], [0],  style='+', label="SC",  fig=fig, ax=ax, color=colors["magenta"])

def NAV_CEL_plot(timeline, DYN_obj, SEN_obj, NAV_obj):
    time_SIM         = timeline["DYN"]["DYN_TIME"]["time_SIM"]
    STRoutflg        = timeline["SEN"]["SEN_STR"]["STRoutflg"]
    SUNdir_SC_mes    = timeline["SEN"]["SEN_STR"]["SUNdir_SC_mes"].T
    EARTHdir_SC_mes  = timeline["SEN"]["SEN_STR"]["EARTHdir_SC_mes"].T
    MOONdir_SC_mes   = timeline["SEN"]["SEN_STR"]["MOONdir_SC_mes"].T
    MARSdir_SC_mes   = timeline["SEN"]["SEN_STR"]["MARSdir_SC_mes"].T
    DEIMOSdir_SC_mes = timeline["SEN"]["SEN_STR"]["DEIMOSdir_SC_mes"].T
    PHOBOSdir_SC_mes = timeline["SEN"]["SEN_STR"]["PHOBOSdir_SC_mes"].T
    STARSdir_SC_mes  = timeline["SEN"]["SEN_STR"]["STARSdir_SC_mes"].T

    NAV_CELoutflg = timeline["NAV"]["NAV_CEL"]["NAV_CELoutflg"]

    SUNangles_mes    = timeline["NAV"]["NAV_CEL"]["SUNangles_mes"]
    EARTHangles_mes  = timeline["NAV"]["NAV_CEL"]["EARTHangles_mes"]
    MOONangles_mes   = timeline["NAV"]["NAV_CEL"]["MOONangles_mes"]
    MARSangles_mes   = timeline["NAV"]["NAV_CEL"]["MARSangles_mes"]
    PHOBOSangles_mes = timeline["NAV"]["NAV_CEL"]["PHOBOSangles_mes"]
    DEIMOSangles_mes = timeline["NAV"]["NAV_CEL"]["DEIMOSangles_mes"]

    SUNangles_mes_deg    = PPC.rad2deg(SUNangles_mes)
    EARTHangles_mes_deg  = PPC.rad2deg(EARTHangles_mes)
    MOONangles_mes_deg   = PPC.rad2deg(MOONangles_mes)
    MARSangles_mes_deg   = PPC.rad2deg(MARSangles_mes)
    PHOBOSangles_mes_deg = PPC.rad2deg(PHOBOSangles_mes)
    DEIMOSangles_mes_deg = PPC.rad2deg(DEIMOSangles_mes)

    SUNsel_STARSdir_SC_mes    = timeline["NAV"]["NAV_CEL"]["SUNsel_STARSdir_SC_mes"].T
    EARTHsel_STARSdir_SC_mes  = timeline["NAV"]["NAV_CEL"]["EARTHsel_STARSdir_SC_mes"].T
    MOONsel_STARSdir_SC_mes   = timeline["NAV"]["NAV_CEL"]["MOONsel_STARSdir_SC_mes"].T
    MARSsel_STARSdir_SC_mes   = timeline["NAV"]["NAV_CEL"]["MARSsel_STARSdir_SC_mes"].T
    PHOBOSsel_STARSdir_SC_mes = timeline["NAV"]["NAV_CEL"]["PHOBOSsel_STARSdir_SC_mes"].T
    DEIMOSsel_STARSdir_SC_mes = timeline["NAV"]["NAV_CEL"]["DEIMOSsel_STARSdir_SC_mes"].T

    # Plot status
    fig, ax = PPC.plot(time_SIM, NAV_CELoutflg, xlabel="time_SIM [s]", ylabel="flag", label="NAV_CELoutflag", title="NAV_CEL output flag")
    PPC.plot(time_SIM, STRoutflg, label="STRoutflag", style='--', fig=fig, ax=ax)

    # Plot angles
    fig, ax = PPC.plot(time_SIM, SUNangles_mes_deg,    xlabel="time_SIM [s]", ylabel="SUNangles_mes [deg]",    label=["$a_1$", "$a_2$", "$a_3$"], title="Line-of-Sight Angle from Sun to Stars")
    fig, ax = PPC.plot(time_SIM, EARTHangles_mes_deg,  xlabel="time_SIM [s]", ylabel="EARTHangles_mes [deg]",  label=["$a_1$", "$a_2$", "$a_3$"], title="Line-of-Sight Angle from Earth to Stars")
    fig, ax = PPC.plot(time_SIM, MOONangles_mes_deg,   xlabel="time_SIM [s]", ylabel="MOONangles_mes [deg]",   label=["$a_1$", "$a_2$", "$a_3$"], title="Line-of-Sight Angle from Moon to Stars")
    fig, ax = PPC.plot(time_SIM, MARSangles_mes_deg,   xlabel="time_SIM [s]", ylabel="MARSangles_mes [deg]",   label=["$a_1$", "$a_2$", "$a_3$"], title="Line-of-Sight Angle from Mars to Stars")
    fig, ax = PPC.plot(time_SIM, DEIMOSangles_mes_deg, xlabel="time_SIM [s]", ylabel="DEIMOSangles_mes [deg]", label=["$a_1$", "$a_2$", "$a_3$"], title="Line-of-Sight Angle from Deimos to Stars")
    fig, ax = PPC.plot(time_SIM, PHOBOSangles_mes_deg, xlabel="time_SIM [s]", ylabel="PHOBOSangles_mes [deg]", label=["$a_1$", "$a_2$", "$a_3$"], title="Line-of-Sight Angle from Phobos to Stars")

    # Plot Selected stars for each body
    fig, ax = PPC.plot(STARSdir_SC_mes[0,:,1], STARSdir_SC_mes[1,:,1], STARSdir_SC_mes[2,:,1], style='.', label="Visible Stars", xlabel="X SC", ylabel="Y SC", zlabel="Z SC", title="Star Field Normalized in SC frame", aspect="equal", color=colors["green"])
    if not PPC.is_nan(SUNangles_mes):
        PPC.plot(SUNsel_STARSdir_SC_mes[0,:,1], SUNsel_STARSdir_SC_mes[1,:,1], SUNsel_STARSdir_SC_mes[2,:,1], style='.', label="Sun selected stars", fig=fig, ax=ax, color=colors["orange"])
        PPC.plot(SUNdir_SC_mes[0,:], SUNdir_SC_mes[1,:], SUNdir_SC_mes[2,:], style='x', label="Sun", fig=fig, ax=ax, color=colors["orange"])
    if not PPC.is_nan(EARTHangles_mes):
        PPC.plot(EARTHsel_STARSdir_SC_mes[0,:,1], EARTHsel_STARSdir_SC_mes[1,:,1], EARTHsel_STARSdir_SC_mes[2,:,1], style='.', label="Earth selected stars",  fig=fig, ax=ax, color=colors["blue"])
        PPC.plot(EARTHdir_SC_mes[0,:], EARTHdir_SC_mes[1,:], EARTHdir_SC_mes[2,:], style='x', label="Earth", fig=fig, ax=ax, color=colors["blue"])
    if not PPC.is_nan(MOONangles_mes):
        PPC.plot(MOONsel_STARSdir_SC_mes[0,:,1], MOONsel_STARSdir_SC_mes[1,:,1], MOONsel_STARSdir_SC_mes[2,:,1], style='.', label="Moon selected stars", fig=fig, ax=ax, color=colors["grey"])
        PPC.plot(MOONdir_SC_mes[0,:], MOONdir_SC_mes[1,:], MOONdir_SC_mes[2,:], style='x', label="Moon", fig=fig, ax=ax, color=colors["grey"])
    if not PPC.is_nan(MARSangles_mes):
        PPC.plot(MARSsel_STARSdir_SC_mes[0,:,1], MARSsel_STARSdir_SC_mes[1,:,1], MARSsel_STARSdir_SC_mes[2,:,1], style='.', label="Mars selected stars", fig=fig, ax=ax, color=colors["red"])
        PPC.plot(MARSdir_SC_mes[0,:], MARSdir_SC_mes[1,:], MARSdir_SC_mes[2,:], style='x', label="Mars", fig=fig, ax=ax, color=colors["red"])
    if not PPC.is_nan(DEIMOSangles_mes):
        PPC.plot(DEIMOSsel_STARSdir_SC_mes[0,:,1], DEIMOSsel_STARSdir_SC_mes[1,:,1], DEIMOSsel_STARSdir_SC_mes[2,:,1], style='.', label="Deimos selected stars", fig=fig, ax=ax, color=colors["lightgrey"])
        PPC.plot(DEIMOSdir_SC_mes[0,:], DEIMOSdir_SC_mes[1,:], DEIMOSdir_SC_mes[2,:], style='x', label="Deimos", fig=fig, ax=ax, color=colors["lightgrey"])
    if not PPC.is_nan(PHOBOSangles_mes):
        PPC.plot(PHOBOSsel_STARSdir_SC_mes[0,:,1], PHOBOSsel_STARSdir_SC_mes[1,:,1], PHOBOSsel_STARSdir_SC_mes[2,:,1], style='.', label="Phobos selected stars", fig=fig, ax=ax, color=colors["darkgrey"])
        PPC.plot(PHOBOSdir_SC_mes[0,:], PHOBOSdir_SC_mes[1,:], PHOBOSdir_SC_mes[2,:], style='x', label="Phobos", fig=fig, ax=ax, color=colors["darkgrey"])

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
    "NAV_CEL"   : NAV_CEL_plot,
}
