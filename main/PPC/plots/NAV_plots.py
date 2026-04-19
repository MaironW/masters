# List all plots to be generated, organized by group name, to be selected on SIM_par

import numpy as np
from PPC import PPC
from Utils import misc
from Utils.constants import CONSTANTS_par

def NAV_EPH_plot(timeline, DYN_obj, SEN_obj, NAV_obj):
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
    fig, ax = PPC.plot(DYN_x, DYN_y, DYN_z, style='.', label="Stars DYN", xlabel="X SSB", ylabel="Y SSB", zlabel="Z SSB", title="Star Field Normalized", aspect="equal", color=PPC.colors["black"])
    PPC.plot(NAV_x, NAV_y, NAV_z, style='.', label="Stars NAV", color=PPC.colors["blue"], fig=fig, ax=ax)

    # 2D sky sphere
    DYN_STARSproj = misc.aitoff_projection(DYN_STARSdir_SSB_reshaped)
    NAV_STARSproj = misc.aitoff_projection(NAV_STARSdir_SSB_reshaped)
    boundary  = misc.aitoff_boundary()
    fig, ax = PPC.plot(DYN_STARSproj[0], DYN_STARSproj[1], style='.', label="DYN Stars", xlabel="Right Ascension [deg]", ylabel="Declination [deg]", title="Star Field - Aitoff Projection", aspect="equal", color=PPC.colors["black"])
    fig, ax = PPC.plot(NAV_STARSproj[0], NAV_STARSproj[1], style='.', label="NAV Stars", color=PPC.colors["blue"], fig=fig, ax=ax)
    PPC.plot(boundary[0], boundary[1], fig=fig, ax=ax)

def NAV_CEL_plot(timeline, DYN_obj, SEN_obj, NAV_obj):
    time_SIM        = timeline["DYN"]["DYN_TIME"]["time_SIM"]
    SCpos_SSB       = timeline["DYN"]["DYN_TRA"]["SCpos_SSB"]
    SCvel_SSB       = timeline["DYN"]["DYN_TRA"]["SCvel_SSB"]
    SEN_STRoutflg   = timeline["SEN"]["SEN_STR"]["SEN_STRoutflg"]
    STARSdir_SC_mes = timeline["SEN"]["SEN_STR"]["STARSdir_SC_mes"]

    NAV_CELoutflg               = timeline["NAV"]["NAV_CEL"]["NAV_CELoutflg"]
    BODYsel_STARdir_SC_mes_list = timeline["NAV"]["NAV_CEL"]["BODYsel_STARdir_SC_mes_list"]

    z = timeline["NAV"]["NAV_CEL"]["z"]
    R = timeline["NAV"]["NAV_CEL"]["R"]

    bodies = [
        dict(name="SUN",    idx=0, color=PPC.colors["orange"]),
        dict(name="MOON",   idx=2, color=PPC.colors["grey"]),
        dict(name="EARTH",  idx=1, color=PPC.colors["blue"]),
        dict(name="DEIMOS", idx=4, color=PPC.colors["darkgrey"]),
        dict(name="PHOBOS", idx=5, color=PPC.colors["lightgrey"]),
        dict(name="MARS",   idx=3, color=PPC.colors["red"]),
    ]

    angles_deg = {}

    for body in bodies:
        idx = body["idx"]
        angles_deg[body["name"]] = (np.cos(z[:, idx])*CONSTANTS_par["rad2deg_cst"])

    active_bodies = [
        body for body in bodies if not misc.is_nan(angles_deg[body["name"]])
    ]

    # Reshape vectors for plot
    STARSdir_SC_mes_reshaped = STARSdir_SC_mes.reshape(-1, 3) # [time * star, direction]

    # Plot status
    fig, ax = PPC.plot(time_SIM, NAV_CELoutflg, xlabel="time_SIM [s]", ylabel="flag", label="NAV_CELoutflag", title="NAV_CEL output flag")
    PPC.plot(time_SIM, SEN_STRoutflg, label="STRoutflag", style='--', fig=fig, ax=ax)

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
    fig, ax = PPC.plot(STARSdir_SC_mes_reshaped[:,0], STARSdir_SC_mes_reshaped[:,1], STARSdir_SC_mes_reshaped[:,2], style='.', label="Visible Stars", xlabel="X SC", ylabel="Y SC", zlabel="Z SC", title="Star Field Normalized in SC frame", aspect="equal", color=PPC.colors["green"])
    PPC.plot([0], [0], [0],  style='+', label="SC",  fig=fig, ax=ax, color=PPC.colors["magenta"])

    for body in active_bodies:
        name = body["name"]
        idx  = body["idx"]
        BODYsel_STARdir_SC_mes = BODYsel_STARdir_SC_mes_list[:,idx,:] # [time, body, direction]
        BODYdir_SC_mes         = timeline["SEN"]["SEN_STR"][f"{name}dir_SC_mes"]

        PPC.plot(BODYsel_STARdir_SC_mes[:,0], BODYsel_STARdir_SC_mes[:,1], BODYsel_STARdir_SC_mes[:,2], style='.', label=f"{name} selected star", fig=fig, ax=ax, color=body["color"])
        PPC.plot(BODYdir_SC_mes[:,0],         BODYdir_SC_mes[:,1],         BODYdir_SC_mes[:,2],         style='x', label=f"{name}",               fig=fig, ax=ax, color=body["color"])

    # Plot selected stars direction over time
    fig, ax1 = PPC.plot(time_SIM, STARSdir_SC_mes[:, :, 0], color=PPC.colors["green"], subplot=(3,1,1), ylabel="x", title=f"Direction of Visible Stars and Reference Selected Star")
    fig, ax2 = PPC.plot(time_SIM, STARSdir_SC_mes[:, :, 1], color=PPC.colors["green"], subplot=(3,1,2), ylabel="y", fig=fig)
    fig, ax3 = PPC.plot(time_SIM, STARSdir_SC_mes[:, :, 2], color=PPC.colors["green"], subplot=(3,1,3), xlabel="time_SIM [s]", ylabel="z", fig=fig)
    for body in active_bodies:
        name = body["name"]
        idx  = body["idx"]
        BODYsel_STARdir_SC_mes = BODYsel_STARdir_SC_mes_list[:,idx,:] # [time, body, direction]
        PPC.plot(time_SIM, BODYsel_STARdir_SC_mes[:, 0], label=f"{name}", color=body["color"], fig=fig, ax=ax1)
        PPC.plot(time_SIM, BODYsel_STARdir_SC_mes[:, 1], label=f"{name}", color=body["color"], fig=fig, ax=ax2)
        PPC.plot(time_SIM, BODYsel_STARdir_SC_mes[:, 2], label=f"{name}", color=body["color"], fig=fig, ax=ax3)

def NAV_PSR_plot(timeline, DYN_obj, SEN_obj, NAV_obj):
    time_SIM      = timeline["DYN"]["DYN_TIME"]["time_SIM"]
    SCpos_SSB     = timeline["DYN"]["DYN_TRA"]["SCpos_SSB"]
    SCvel_SSB     = timeline["DYN"]["DYN_TRA"]["SCvel_SSB"]
    SEN_PSRoutflg = timeline["SEN"]["SEN_PSR"]["SEN_PSRoutflg"]
    NAV_PSRoutflg = timeline["NAV"]["NAV_PSR"]["NAV_PSRoutflg"]

    n_pulsars  = NAV_obj.NAV_PSR.par["n_pulsars"]
    PULSARname = NAV_obj.NAV_PSR.par["name"]

    z = timeline["NAV"]["NAV_PSR"]["z"]
    R = timeline["NAV"]["NAV_PSR"]["R"]

    # Plot status
    fig, ax = PPC.plot(time_SIM, NAV_PSRoutflg, xlabel="time_SIM [s]", ylabel="flag", label="NAV_PSRoutflag", title="NAV_PSR output flag")
    PPC.plot(time_SIM, SEN_PSRoutflg, label="PSRoutflag", style='--', fig=fig, ax=ax)

    # Plot times and covariance
    R_diag  = np.diagonal(R, axis1=1, axis2=2)
    sigma_z = np.sqrt(R_diag)

    fig, ax1 = PPC.plot([], [], ylabel="OBTdt_TDB_mes [s]", title="TOA delay on SC", subplot=(2,1,1))
    fig, ax2 = PPC.plot([], [], xlabel="time_SIM [s]", title="Standard deviation", ylabel="Covariance [s]", fig=fig, subplot=(2,1,2))
    for i in range(n_pulsars):
        PPC.plot(time_SIM, z[:, i],       label=f"{PULSARname[i]} meas.", fig=fig, ax=ax1)
        PPC.plot(time_SIM, sigma_z[:, i], label=f"{PULSARname[i]} σ",     fig=fig, ax=ax2)

    # Plot predicted measurement, assuming the true spacecraft position as the state + innovation
    x_true = np.hstack((SCpos_SSB, SCvel_SSB))
    n_iter = len(x_true)
    h_hist = np.full((n_iter, n_pulsars), np.nan)
    H_hist = np.zeros((n_iter, n_pulsars, 6))
    # Iterate from idx=1 onwards, because we need to use x_true[k-1] to compute h(x)
    for k in range(n_iter):
        if NAV_PSRoutflg[k] == 1:
            # First update NAV_PSR state, otherwise h(x) will be computed for the last (already computed) state
            NAV_obj.NAV_PSR.state["SSBpos_SUN_ref"] = timeline["NAV"]["NAV_PSR"]["SSBpos_SUN_ref"][k]
            h_hist[k] = NAV_obj.NAV_PSR.h(x_true[k])
            H_hist[k] = NAV_obj.NAV_PSR.H(x_true[k])
    fig, ax1 = PPC.plot([], [], ylabel="h(x)", title="Predicted Measurement h(x)", subplot=(2,1,1))
    fig, ax2 = PPC.plot([], [], xlabel="time_SIM [s]", ylabel="z - h(x)", title="Innovation z - h(x)", fig=fig, subplot=(2,1,2))
    innov = z - h_hist
    for i in range(n_pulsars):
        name = PULSARname[i]
        PPC.plot(time_SIM, h_hist[:, i], label=f"h(x) {name}", fig=fig, ax=ax1)
        PPC.plot(time_SIM, z[:, i], label=f"z {name}", style='--', fig=fig, ax=ax1)
        PPC.plot(time_SIM, innov[:, i], label=f"{name}", xlabel="time_SIM [s]", fig=fig, ax=ax2)

def NAV_CMB_plot(timeline, DYN_obj, SEN_obj, NAV_obj):
    time_SIM      = timeline["DYN"]["DYN_TIME"]["time_SIM"]
    SCpos_SSB     = timeline["DYN"]["DYN_TRA"]["SCpos_SSB"]
    SCvel_SSB     = timeline["DYN"]["DYN_TRA"]["SCvel_SSB"]
    SEN_STRoutflg = timeline["SEN"]["SEN_STR"]["SEN_STRoutflg"]
    SEN_CMBoutflg = timeline["SEN"]["SEN_CMB"]["SEN_CMBoutflg"]
    NAV_CMBoutflg = timeline["NAV"]["NAV_CMB"]["NAV_CMBoutflg"]

    z = timeline["NAV"]["NAV_CMB"]["z"]
    R = timeline["NAV"]["NAV_CMB"]["R"]

    # Plot status
    fig, ax = PPC.plot(time_SIM, NAV_CMBoutflg, xlabel="time_SIM [s]", ylabel="flag", label="NAV_CMBoutflag", title="NAV_CMB output flag")
    PPC.plot(time_SIM, SEN_CMBoutflg, label="CMBoutflag", style='--', fig=fig, ax=ax)
    PPC.plot(time_SIM, SEN_STRoutflg, label="STRoutflag", style='--', fig=fig, ax=ax)

    # Plot times and covariance
    R_diag  = np.diagonal(R, axis1=1, axis2=2)
    sigma_z = np.sqrt(R_diag)

    fig, ax1 = PPC.plot([], [], ylabel="CMBR Temperature Dipole [K]", title="CMBR Temperature Dipole", subplot=(2,1,1))
    fig, ax2 = PPC.plot([], [], xlabel="time_SIM [s]", title="Standard deviation", ylabel="Covariance [K]", fig=fig, subplot=(2,1,2))
    for i in range(3):
        PPC.plot(time_SIM, z[:, i],       label=f"CMB{i} meas.", fig=fig, ax=ax1)
        PPC.plot(time_SIM, sigma_z[:, i], label=f"CMB{i} σ",     fig=fig, ax=ax2)

    # Plot predicted measurement, assuming the true spacecraft velocity as the state + innovation
    x_true = np.hstack((SCpos_SSB, SCvel_SSB))
    n_iter = len(x_true)
    h_hist = np.full((n_iter, 3), np.nan)
    H_hist = np.zeros((n_iter, 3, 6))
    # Iterate from idx=1 onwards, because we need to use x_true[k-1] to compute h(x)
    for k in range(n_iter):
        if NAV_CMBoutflg[k] == 1:
            # First update NAV_CMB state, otherwise h(x) will be computed for the last (already computed) state
            NAV_obj.NAV_CMB.state["BOFq_SSB_ref"] = timeline["NAV"]["NAV_CMB"]["BOFq_SSB_ref"][k]
            h_hist[k] = NAV_obj.NAV_CMB.h(x_true[k])
            H_hist[k] = NAV_obj.NAV_CMB.H(x_true[k])
    fig, ax1 = PPC.plot([], [], ylabel="h(x)", title="Predicted Measurement h(x)", subplot=(2,1,1))
    fig, ax2 = PPC.plot([], [], xlabel="time_SIM [s]", ylabel="z - h(x)", title="Innovation z - h(x)", fig=fig, subplot=(2,1,2))
    innov = z - h_hist
    for i in range(3):
        PPC.plot(time_SIM, h_hist[:, i], label=f"h(x) CMB{i}", fig=fig, ax=ax1)
        PPC.plot(time_SIM, z[:, i], label=f"z CMB{i}", style='--', fig=fig, ax=ax1)
        PPC.plot(time_SIM, innov[:, i], label=f"CMB{i}", xlabel="time_SIM [s]", fig=fig, ax=ax2)

def NAV_EKF_plot(timeline, DYN_obj, SEN_obj, NAV_obj):
    time_SIM      = timeline["DYN"]["DYN_TIME"]["time_SIM"]
    SCpos_SSB     = timeline["DYN"]["DYN_TRA"]["SCpos_SSB"]
    SCvel_SSB     = timeline["DYN"]["DYN_TRA"]["SCvel_SSB"]
    SEN_STRoutflg = timeline["SEN"]["SEN_STR"]["SEN_STRoutflg"]
    SEN_CMBoutflg = timeline["SEN"]["SEN_CMB"]["SEN_CMBoutflg"]
    NAV_CMBoutflg = timeline["NAV"]["NAV_CMB"]["NAV_CMBoutflg"]

    x_est = timeline["NAV"]["NAV_EKF"]["x_est"]
    SCpos_SSB_est = x_est[:, 0:3]
    SCvel_SSB_est = x_est[:, 3:6]

    fig, ax1 = PPC.plot([], [], ylabel="SCpos_SSB [km]", title="Estimated position", subplot=(2,1,1))
    fig, ax2 = PPC.plot([], [], xlabel="time_SIM [s]", ylabel="SCvel_SSB [km/s]", title="Estimated velocity", fig=fig, subplot=(2,1,2))
    PPC.plot(time_SIM, SCpos_SSB,     label=f"SCpos_SSB",     fig=fig, ax=ax1)
    PPC.plot(time_SIM, SCpos_SSB_est, label=f"SCpos_SSB_est", fig=fig, ax=ax1)
    PPC.plot(time_SIM, SCvel_SSB,     label=f"SCvel_SSB",     fig=fig, ax=ax2)
    PPC.plot(time_SIM, SCvel_SSB_est, label=f"SCvel_SSB_est", fig=fig, ax=ax2)

    fig, ax1 = PPC.plot([], [], ylabel="SCpos_SSB [km]", title="Estimated position error", subplot=(2,1,1))
    fig, ax2 = PPC.plot([], [], xlabel="time_SIM [s]", ylabel="SCvel_SSB [km/s]", title="Estimated velocity error", fig=fig, subplot=(2,1,2))
    PPC.plot(time_SIM, SCpos_SSB - SCpos_SSB_est, label=f"SCpos_SSB", fig=fig, ax=ax1)
    PPC.plot(time_SIM, SCvel_SSB - SCvel_SSB_est, label=f"SCvel_SSB", fig=fig, ax=ax2)

NAV_plots = {
    "NAV_EPH" : NAV_EPH_plot,
    "NAV_STR" : NAV_STR_plot,
    "NAV_CEL" : NAV_CEL_plot,
    "NAV_PSR" : NAV_PSR_plot,
    "NAV_CMB" : NAV_CMB_plot,
    "NAV_EKF" : NAV_EKF_plot,
}
