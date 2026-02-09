# List all plots to be generated, organized by group name, to be selected on SIM_par

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

def SEN_STR_plot(timeline, DYN_obj, SEN_obj, NAV_obj):
    time_SIM     = timeline["DYN"]["DYN_TIME"]["time_SIM"]
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

    # Compare STR and SIM times
    fig, ax = PPC.plot(time_SIM, time_SIM, xlabel="time_SIM [s]", ylabel="time [s]", label="time_SIM", title="STR Time")
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


def NAV_CEL_plot(timeline, DYN_obj, SEN_obj, NAV_obj):
    time_SIM         = timeline["DYN"]["DYN_TIME"]["time_SIM"]
    STRoutflg        = timeline["SEN"]["SEN_STR"]["STRoutflg"]
    SUNdir_SC_mes    = timeline["SEN"]["SEN_STR"]["SUNdir_SC_mes"]
    EARTHdir_SC_mes  = timeline["SEN"]["SEN_STR"]["EARTHdir_SC_mes"]
    MOONdir_SC_mes   = timeline["SEN"]["SEN_STR"]["MOONdir_SC_mes"]
    MARSdir_SC_mes   = timeline["SEN"]["SEN_STR"]["MARSdir_SC_mes"]
    DEIMOSdir_SC_mes = timeline["SEN"]["SEN_STR"]["DEIMOSdir_SC_mes"]
    PHOBOSdir_SC_mes = timeline["SEN"]["SEN_STR"]["PHOBOSdir_SC_mes"]
    STARSdir_SC_mes  = timeline["SEN"]["SEN_STR"]["STARSdir_SC_mes"]

    NAV_CELoutflg = timeline["NAV"]["NAV_CEL"]["NAV_CELoutflg"]

    SUNangles_mes    = timeline["NAV"]["NAV_CEL"]["SUNangles_mes"]
    EARTHangles_mes  = timeline["NAV"]["NAV_CEL"]["EARTHangles_mes"]
    MOONangles_mes   = timeline["NAV"]["NAV_CEL"]["MOONangles_mes"]
    MARSangles_mes   = timeline["NAV"]["NAV_CEL"]["MARSangles_mes"]
    PHOBOSangles_mes = timeline["NAV"]["NAV_CEL"]["PHOBOSangles_mes"]
    DEIMOSangles_mes = timeline["NAV"]["NAV_CEL"]["DEIMOSangles_mes"]

    SUNangles_mes_deg    = SUNangles_mes*CONSTANTS_par["rad2deg_cst"]
    EARTHangles_mes_deg  = EARTHangles_mes*CONSTANTS_par["rad2deg_cst"]
    MOONangles_mes_deg   = MOONangles_mes*CONSTANTS_par["rad2deg_cst"]
    MARSangles_mes_deg   = MARSangles_mes*CONSTANTS_par["rad2deg_cst"]
    PHOBOSangles_mes_deg = PHOBOSangles_mes*CONSTANTS_par["rad2deg_cst"]
    DEIMOSangles_mes_deg = DEIMOSangles_mes*CONSTANTS_par["rad2deg_cst"]

    SUNsel_STARSdir_SC_mes    = timeline["NAV"]["NAV_CEL"]["SUNsel_STARSdir_SC_mes"]
    EARTHsel_STARSdir_SC_mes  = timeline["NAV"]["NAV_CEL"]["EARTHsel_STARSdir_SC_mes"]
    MOONsel_STARSdir_SC_mes   = timeline["NAV"]["NAV_CEL"]["MOONsel_STARSdir_SC_mes"]
    MARSsel_STARSdir_SC_mes   = timeline["NAV"]["NAV_CEL"]["MARSsel_STARSdir_SC_mes"]
    PHOBOSsel_STARSdir_SC_mes = timeline["NAV"]["NAV_CEL"]["PHOBOSsel_STARSdir_SC_mes"]
    DEIMOSsel_STARSdir_SC_mes = timeline["NAV"]["NAV_CEL"]["DEIMOSsel_STARSdir_SC_mes"]

    # Reshape vectors for plot
    STARSdir_SC_mes_reshaped           = STARSdir_SC_mes.reshape(-1, 3)           # [time * star, direction]
    SUNsel_STARSdir_SC_mes_reshaped    = SUNsel_STARSdir_SC_mes.reshape(-1, 3)    # [time * star, direction]
    EARTHsel_STARSdir_SC_mes_reshaped  = EARTHsel_STARSdir_SC_mes.reshape(-1, 3)  # [time * star, direction]
    MOONsel_STARSdir_SC_mes_reshaped   = MOONsel_STARSdir_SC_mes.reshape(-1, 3)   # [time * star, direction]
    MARSsel_STARSdir_SC_mes_reshaped   = MARSsel_STARSdir_SC_mes.reshape(-1, 3)   # [time * star, direction]
    PHOBOSsel_STARSdir_SC_mes_reshaped = PHOBOSsel_STARSdir_SC_mes.reshape(-1, 3) # [time * star, direction]
    DEIMOSsel_STARSdir_SC_mes_reshaped = DEIMOSsel_STARSdir_SC_mes.reshape(-1, 3) # [time * star, direction]

    # Plot status
    fig, ax = PPC.plot(time_SIM, NAV_CELoutflg, xlabel="time_SIM [s]", ylabel="flag", label="NAV_CELoutflag", title="NAV_CEL output flag")
    PPC.plot(time_SIM, STRoutflg, label="STRoutflag", style='--', fig=fig, ax=ax)

    # Plot angles
    if not PPC.is_nan(SUNangles_mes_deg):    PPC.plot(time_SIM, SUNangles_mes_deg,    xlabel="time_SIM [s]", ylabel="SUNangles_mes [deg]",    label=["$a_1$", "$a_2$", "$a_3$"], title="Line-of-Sight Angle from Sun to Stars")
    if not PPC.is_nan(MOONangles_mes_deg):   PPC.plot(time_SIM, MOONangles_mes_deg,   xlabel="time_SIM [s]", ylabel="MOONangles_mes [deg]",   label=["$a_1$", "$a_2$", "$a_3$"], title="Line-of-Sight Angle from Moon to Stars")
    if not PPC.is_nan(EARTHangles_mes_deg):  PPC.plot(time_SIM, EARTHangles_mes_deg,  xlabel="time_SIM [s]", ylabel="EARTHangles_mes [deg]",  label=["$a_1$", "$a_2$", "$a_3$"], title="Line-of-Sight Angle from Earth to Stars")
    if not PPC.is_nan(DEIMOSangles_mes_deg): PPC.plot(time_SIM, DEIMOSangles_mes_deg, xlabel="time_SIM [s]", ylabel="DEIMOSangles_mes [deg]", label=["$a_1$", "$a_2$", "$a_3$"], title="Line-of-Sight Angle from Deimos to Stars")
    if not PPC.is_nan(PHOBOSangles_mes_deg): PPC.plot(time_SIM, PHOBOSangles_mes_deg, xlabel="time_SIM [s]", ylabel="PHOBOSangles_mes [deg]", label=["$a_1$", "$a_2$", "$a_3$"], title="Line-of-Sight Angle from Phobos to Stars")
    if not PPC.is_nan(MARSangles_mes_deg):   PPC.plot(time_SIM, MARSangles_mes_deg,   xlabel="time_SIM [s]", ylabel="MARSangles_mes [deg]",   label=["$a_1$", "$a_2$", "$a_3$"], title="Line-of-Sight Angle from Mars to Stars")

    # Plot Selected stars for each body
    fig, ax = PPC.plot(STARSdir_SC_mes_reshaped[:,0], STARSdir_SC_mes_reshaped[:,1], STARSdir_SC_mes_reshaped[:,2], style='.', label="Visible Stars", xlabel="X SC", ylabel="Y SC", zlabel="Z SC", title="Star Field Normalized in SC frame", aspect="equal", color=colors["green"])
    if not PPC.is_nan(SUNangles_mes):
        PPC.plot(SUNsel_STARSdir_SC_mes_reshaped[:,0], SUNsel_STARSdir_SC_mes_reshaped[:,1], SUNsel_STARSdir_SC_mes_reshaped[:,2], style='.', label="Sun selected stars", fig=fig, ax=ax, color=colors["orange"])
        PPC.plot(SUNdir_SC_mes[:,0],                   SUNdir_SC_mes[:,1],                   SUNdir_SC_mes[:,2],                   style='x', label="Sun",                fig=fig, ax=ax, color=colors["orange"])
    if not PPC.is_nan(MOONangles_mes):
        PPC.plot(MOONsel_STARSdir_SC_mes_reshaped[:,0], MOONsel_STARSdir_SC_mes_reshaped[:,1], MOONsel_STARSdir_SC_mes_reshaped[:,2], style='.', label="Moon selected stars", fig=fig, ax=ax, color=colors["grey"])
        PPC.plot(MOONdir_SC_mes[:,0],                   MOONdir_SC_mes[:,1],                   MOONdir_SC_mes[:,2],                   style='x', label="Moon",                fig=fig, ax=ax, color=colors["grey"])
    if not PPC.is_nan(EARTHangles_mes):
        PPC.plot(EARTHsel_STARSdir_SC_mes_reshaped[:,0], EARTHsel_STARSdir_SC_mes_reshaped[:,1], EARTHsel_STARSdir_SC_mes_reshaped[:,2], style='.', label="Earth selected stars",  fig=fig, ax=ax, color=colors["blue"])
        PPC.plot(EARTHdir_SC_mes[:,0],                   EARTHdir_SC_mes[:,1],                   EARTHdir_SC_mes[:,2],                   style='x', label="Earth",                 fig=fig, ax=ax, color=colors["blue"])
    if not PPC.is_nan(DEIMOSangles_mes):
        PPC.plot(DEIMOSsel_STARSdir_SC_mes_reshaped[:,0], DEIMOSsel_STARSdir_SC_mes_reshaped[:,1], DEIMOSsel_STARSdir_SC_mes_reshaped[:,2], style='.', label="Deimos selected stars", fig=fig, ax=ax, color=colors["lightgrey"])
        PPC.plot(DEIMOSdir_SC_mes[:,0],                   DEIMOSdir_SC_mes[:,1],                   DEIMOSdir_SC_mes[:,2],                   style='x', label="Deimos",                fig=fig, ax=ax, color=colors["lightgrey"])
    if not PPC.is_nan(PHOBOSangles_mes):
        PPC.plot(PHOBOSsel_STARSdir_SC_mes_reshaped[:,0], PHOBOSsel_STARSdir_SC_mes_reshaped[:,1], PHOBOSsel_STARSdir_SC_mes_reshaped[:,2], style='.', label="Phobos selected stars", fig=fig, ax=ax, color=colors["darkgrey"])
        PPC.plot(PHOBOSdir_SC_mes[:,0],                   PHOBOSdir_SC_mes[:,1],                   PHOBOSdir_SC_mes[:,2],                   style='x', label="Phobos",                fig=fig, ax=ax, color=colors["darkgrey"])
    if not PPC.is_nan(MARSangles_mes):
        PPC.plot(MARSsel_STARSdir_SC_mes_reshaped[:,0], MARSsel_STARSdir_SC_mes_reshaped[:,1], MARSsel_STARSdir_SC_mes_reshaped[:,2], style='.', label="Mars selected stars", fig=fig, ax=ax, color=colors["red"])
        PPC.plot(MARSdir_SC_mes[:,0],                   MARSdir_SC_mes[:,1],                   MARSdir_SC_mes[:,2],                   style='x', label="Mars",                fig=fig, ax=ax, color=colors["red"])
    PPC.plot([0], [0], [0],  style='+', label="SC",  fig=fig, ax=ax, color=colors["magenta"])

    # Plot selected stars direction over time
    fig, ax1 = PPC.plot(time_SIM, STARSdir_SC_mes[:, :, 0], color=colors["green"], subplot=(3,1,1), xlabel="time_SIM [s]", ylabel="x", title="Direction of Visible Stars and Selected Stars for Mars CeleNav")
    PPC.plot(time_SIM, MARSsel_STARSdir_SC_mes[:, :, 0], color=colors["red"], fig=fig, ax=ax1)

    fig, ax2 = PPC.plot(time_SIM, STARSdir_SC_mes[:, :, 1], color=colors["green"], subplot=(3,1,2), fig=fig, ylabel='y')
    PPC.plot(time_SIM, MARSsel_STARSdir_SC_mes[:, :, 1], color=colors["red"], fig=fig, ax=ax2)

    fig, ax3 = PPC.plot(time_SIM, STARSdir_SC_mes[:, :, 2], color=colors["green"], subplot=(3,1,3), fig=fig, ylabel='z')
    PPC.plot(time_SIM, MARSsel_STARSdir_SC_mes[:, :, 2], color=colors["red"], fig=fig, ax=ax3)

def DYN_PSR_plot(timeline, DYN_obj, SEN_obj, NAV_obj):
    STARSdir_SSB   = timeline["DYN"]["DYN_STR"]["STARSdir_SSB"][0] # [time, star, direction]
    PULSARSdir_SSB = timeline["DYN"]["DYN_PSR"]["PULSARSdir_SSB"][0] # [time, pulsar, direction]

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
    boundary  = PPC.aitoff_boundary()
    fig, ax = PPC.plot(STARSproj[0], STARSproj[1], style='.', label="Stars", xlabel="Right Ascension [deg]", ylabel="Declination [deg]", title="Star Field - Aitoff Projection", aspect="equal", color=colors["black"])
    PPC.plot(PULSARSproj[0], PULSARSproj[1], style='x', label="Pulsars", color=colors["purple"], fig=fig, ax=ax)
    PPC.plot(boundary[0], boundary[1], fig=fig, ax=ax)

PPC_plots = {
    "DYN_TIME"  : DYN_TIME_plot,
    "DYN_SUN"   : DYN_SUN_plot,
    "DYN_EARTH" : DYN_EARTH_plot,
    "DYN_MARS"  : DYN_MARS_plot,
    "DYN_GRV"   : DYN_GRV_plot,
    "DYN_TRA"   : DYN_TRA_plot,
    "DYN_ATT"   : DYN_ATT_plot,
    "DYN_STR"   : DYN_STR_plot,
    "DYN_PSR"   : DYN_PSR_plot,
    "SEN_STR"   : SEN_STR_plot,
    "NAV_CEL"   : NAV_CEL_plot,
}
