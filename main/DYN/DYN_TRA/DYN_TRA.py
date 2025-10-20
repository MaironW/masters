# Level 2 Module DYN_TRA
# Simulates the propagation of the spacecraft orbital states for the simulation

from DYN.DYN_GRV import DYN_GRV
from .DYN_TRA_par import DYN_TRA_par
from Utils.constants import CONSTANTS_par
from Utils import quaternions

# Module output dictionary
DYN_TRA_out = {
    "SCpos_TER" : DYN_TRA_par["SCpos_TER_ini"],
    "SCpos_ECI" : DYN_TRA_par["SCpos_ECI_ini"],
    "SCpos_MAR" : DYN_TRA_par["SCpos_MAR_ini"],
    "SCpos_SUN" : DYN_TRA_par["SCpos_SUN_ini"],
    "SCpos_SSB" : DYN_TRA_par["SCpos_SSB_ini"],
    "SCvel_ECI" : DYN_TRA_par["SCvel_ECI_ini"],
    "SCvel_SSB" : DYN_TRA_par["SCvel_SSB_ini"]
}

# Module main function
def run(state_prev, DYN_TIME_out, DYN_EARTH_out):
    dt = DYN_TIME_out["dt"]

    # Integrate coupled equations
    state_next = rk4_step(state_prev, dt, derivatives, DYN_EARTH_out)

    TERq_ECI = DYN_EARTH_out["TERq_ECI"]
    SCpos_ECI = state_next["SCpos_ECI"]
    SCvel_ECI = state_next["SCvel_ECI"]
    SCpos_TER = quaternions.qvecrot(SCpos_ECI, TERq_ECI)

    SCpos_SSB = SCpos_ECI + DYN_EARTH_out["EARTHpos_SSB"]
    SCvel_SSB = SCvel_ECI + DYN_EARTH_out["EARTHvel_SSB"]

    DYN_TRA_out["SCpos_TER"] = SCpos_TER
    DYN_TRA_out["SCpos_ECI"] = SCpos_ECI
    DYN_TRA_out["SCvel_ECI"] = SCvel_ECI
    DYN_TRA_out["SCpos_SSB"] = SCpos_SSB
    DYN_TRA_out["SCvel_SSB"] = SCvel_SSB

    return dict(DYN_TRA_out)

# Module computation of derivatives to be integrated
def derivatives(state, DYN_EARTH_out):
    SCpos_ECI = state["SCpos_ECI"]
    SCvel_ECI = state["SCvel_ECI"]

    TERq_ECI = DYN_EARTH_out["TERq_ECI"]
    SCpos_TER = quaternions.qvecrot(SCpos_ECI, TERq_ECI)

    # Build a minimal DYN_TRA_out to feed into DYN_GRV
    DYN_TRA_out_tmp = {
        "SCpos_TER" : SCpos_TER,
    }

    # Compute gravitational acceleration
    DYN_GRV_out = DYN_GRV.run(DYN_EARTH_out, DYN_TRA_out_tmp)
    grvacc_ECI = DYN_GRV_out["grvacc_ECI"]

    # Return derivatives
    dSCpos_ECI = SCvel_ECI
    dSCvel_ECI = grvacc_ECI

    return {
        "dSCpos_ECI" : dSCpos_ECI,
        "dSCvel_ECI" : dSCvel_ECI
    }

# Runge-Kutta 4 integrator.
# TODO: Set this as a generic function in utils to be used on other modules
def rk4_step(state, dt, derivatives, *args):
    k1 = derivatives(state, *args)

    state_k2 = {
        "SCpos_ECI" : state["SCpos_ECI"] + 0.5*dt*k1["dSCpos_ECI"],
        "SCvel_ECI" : state["SCvel_ECI"] + 0.5*dt*k1["dSCvel_ECI"]
    }
    k2 = derivatives(state_k2, *args)

    state_k3 = {
        "SCpos_ECI" : state["SCpos_ECI"] + 0.5*dt*k2["dSCpos_ECI"],
        "SCvel_ECI" : state["SCvel_ECI"] + 0.5*dt*k2["dSCvel_ECI"]
    }
    k3 = derivatives(state_k3, *args)

    state_k4 = {
        "SCpos_ECI" : state["SCpos_ECI"] + dt*k3["dSCpos_ECI"],
        "SCvel_ECI" : state["SCvel_ECI"] + dt*k3["dSCvel_ECI"]
    }
    k4 = derivatives(state_k4, *args)

    SCpos_ECI = state["SCpos_ECI"] + dt/6.0 * (k1["dSCpos_ECI"] + 2*k2["dSCpos_ECI"] + 2*k3["dSCpos_ECI"] + k4["dSCpos_ECI"])
    SCvel_ECI = state["SCvel_ECI"] + dt/6.0 * (k1["dSCvel_ECI"] + 2*k2["dSCvel_ECI"] + 2*k3["dSCvel_ECI"] + k4["dSCvel_ECI"])

    return {
        "SCpos_ECI" : SCpos_ECI,
        "SCvel_ECI" : SCvel_ECI,
    }
