
# List all plots to be generated for DYN modules, organized by group name, to be selected on SIM_par

from PPC import PPC

def DYN_TIME_plot(timeline, DYN_obj, SEN_obj, NAV_obj):
    PPC.plot(timeline["DYN"]["DYN_TIME"]["time_SIM"], timeline["DYN"]["DYN_TIME"]["time_TDB"], xlabel="time_SIM [s]", ylabel="time_TDB [s]", title="Time")

def DYN_EPH_plot(timeline, DYN_obj, SEN_obj, NAV_obj):
    time_SIM = timeline["DYN"]["DYN_TIME"]["time_SIM"]

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
        # Plot reference position from DYN
        fig, ax1 = PPC.plot(time_SIM, timeline["DYN"]["DYN_EPH"][f"{name}pos_SSB"], ylabel=f"{name}pos_SSB [km]", label=["x","y","z"], title=f"{name}pos_SSB", subplot=(2,1,1))
        # Plot reference velocity from DYN
        fig, ax2 = PPC.plot(time_SIM, timeline["DYN"]["DYN_EPH"][f"{name}vel_SSB"], ylabel=f"{name}vel_SSB [km]", label=["x","y","z"], title=f"{name}vel_SSB", fig=fig, subplot=(2,1,2))

def DYN_GRV_plot(timeline, DYN_obj, SEN_obj, NAV_obj):
    # Gravity acceleration on SSB frame
    PPC.plot(timeline["DYN"]["DYN_TIME"]["time_SIM"], timeline["DYN"]["DYN_GRV"]["grvacc_SSB"], ylabel="grvacc_SSB [km/s^2]", label=["x","y","z"], title="grvacc_SSB")

def DYN_TRA_plot(timeline, DYN_obj, SEN_obj, NAV_obj):
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
    PPC.plot(timeline["DYN"]["DYN_TRA"]["SCpos_SSB"][:,0], timeline["DYN"]["DYN_TRA"]["SCpos_SSB"][:,1], timeline["DYN"]["DYN_TRA"]["SCpos_SSB"][:,2], label="SC", fig=fig, ax=ax, color=PPC.colors["magenta"])

    # 2D SC and Moon trajectories around Earth
    MOONpos_ECI = timeline["DYN"]["DYN_EPH"]["EARTHpos_SSB"] - timeline["DYN"]["DYN_EPH"]["MOONpos_SSB"]
    fig, ax = PPC.plot([0], [0], style='o', label="Earth", xlabel="X ECI [km]", ylabel="Y ECI [km]", title="Trajectories around Earth", aspect="equal", color=PPC.colors["blue"])
    PPC.plot(timeline["DYN"]["DYN_TRA"]["SCpos_ECI"][:,0], timeline["DYN"]["DYN_TRA"]["SCpos_ECI"][:,1], label="SC",   fig=fig, ax=ax, color=PPC.colors["magenta"])
    PPC.plot(MOONpos_ECI[:,0],                             MOONpos_ECI[:,1],                             label="Moon", fig=fig, ax=ax, color=PPC.colors["grey"])

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
    PPC.plot(x, y, z, style='.', label="Stars", xlabel="X SSB", ylabel="Y SSB", zlabel="Z SSB", title="Star Field Normalized", aspect="equal", color=PPC.colors["black"])

    # 2D sky sphere
    STARSproj = PPC.aitoff_projection(STARSdir_SSB_reshaped)
    boundary  = PPC.aitoff_boundary()
    fig, ax = PPC.plot(STARSproj[0], STARSproj[1], style='.', label="Stars", xlabel="Right Ascension [deg]", ylabel="Declination [deg]", title="Star Field - Aitoff Projection", aspect="equal", color=PPC.colors["black"])
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
    fig, ax = PPC.plot(x_star, y_star, z_star, style='.', label="Stars",   xlabel="X SSB", ylabel="Y SSB", zlabel="Z SSB", title="Star Field Normalized", aspect="equal", color=PPC.colors["black"])
    PPC.plot(x_pulsar, y_pulsar, z_pulsar, style='x', label="Pulsars", color=PPC.colors["purple"], fig=fig, ax=ax)

    # 2D sky sphere
    STARSproj   = PPC.aitoff_projection(STARSdir_SSB_reshaped)
    PULSARSproj = PPC.aitoff_projection(PULSARSdir_SSB_reshaped)
    boundary    = PPC.aitoff_boundary()
    fig, ax = PPC.plot(STARSproj[0], STARSproj[1], style='.', label="Stars", xlabel="Right Ascension [deg]", ylabel="Declination [deg]", title="Star Field - Aitoff Projection", aspect="equal", color=PPC.colors["black"])
    PPC.plot(PULSARSproj[0], PULSARSproj[1], style='x', label="Pulsars", color=PPC.colors["purple"], fig=fig, ax=ax)
    PPC.plot(boundary[0], boundary[1], fig=fig, ax=ax)

    # Delays for TOAs between SC and SSB
    fig, ax1 = PPC.plot([], [], ylabel="SCdt_SSB [s]", title="True TOA delay on SC", subplot=(3,1,1))
    fig, ax2 = PPC.plot([], [], ylabel="roemer_delay [s]", fig=fig, subplot=(3,1,2))
    fig, ax3 = PPC.plot([], [], xlabel="time_SIM [s]", ylabel="shapiro_delay [s]", fig=fig, subplot=(3,1,3))
    for i in range(n_pulsars):
        PPC.plot(time_SIM, SCdt_SSB[:,i], label=PULSARname[i], xlabel="time_SIM [s]", ylabel="SCdt_SSB [s]", title="True TOA delay on SC", fig=fig, ax=ax1)
        PPC.plot(time_SIM, roemer_delay[:,i], label=PULSARname[i], ylabel="roemer_delay [s]", fig=fig, ax=ax2)
        PPC.plot(time_SIM, shapiro_delay[:,i], label=PULSARname[i], ylabel="shapiro_delay [s]", fig=fig, ax=ax3)

DYN_plots = {
    "DYN_TIME" : DYN_TIME_plot,
    "DYN_EPH"  : DYN_EPH_plot,
    "DYN_GRV"  : DYN_GRV_plot,
    "DYN_TRA"  : DYN_TRA_plot,
    "DYN_ATT"  : DYN_ATT_plot,
    "DYN_STR"  : DYN_STR_plot,
    "DYN_PSR"  : DYN_PSR_plot,
}
