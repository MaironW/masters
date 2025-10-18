# Level 2 Module DYN_TRA
# Simulates the propagation of the spacecraft orbital states for the simulation

from .DYN_TRA_par import DYN_TRA_par
from Utils.constants import CONSTANTS_par

# Module output dictionary
DYN_TRA_out = {
    "SCpos_TER" : DYN_TRA_par["SCpos_TER_ini"],
    "SCpos_MAR" : DYN_TRA_par["SCpos_MAR_ini"],
    "SCpos_SUN" : DYN_TRA_par["SCpos_SUN_ini"],
    "SCpos_SSB" : DYN_TRA_par["SCpos_SSB_ini"],
    "SCvel_SSB" : DYN_TRA_par["SCvel_SSB_ini"]
}

# Module main function
def run(DYN_TIME_out, DYN_SUN_out, DYN_EARTH_out, DYN_MARS_out, DYN_ATT_out, DYN_GRV_out):
    # Get gravity expressed in each body inertial frame
    grvacc_ECI = DYN_GRV_out["grvacc_ECI"] # [km/s^2]
    grvacc_MCI = DYN_GRV_out["grvacc_MCI"] # [km/s^2]
    grvacc_SSB = DYN_GRV_out["grvacc_SSB"] # [km/s^2]

    # Compute the acceleration applied on the spacecraft (gravity + environmental forces(TBD))
    SCacc_ECI = grvacc_ECI # [km/s^2]
    SCacc_MCI = grvacc_MCI # [km/s^2]
    SCacc_SSB = grvacc_SSB # [km/s^2]

    # Integrate to get position and velocity
    dt = DYN_TIME_out["dt"]
    SCpos_ECI, SCvel_ECI = rk4_step(dt, SCpos_ECI, SCvel_ECI, SCacc_ECI)

    DYN_TRA_out["SCpos_TER"] = DYN_TRA_par["SCpos_TER_ini"]
    DYN_TRA_out["SCpos_MAR"] = DYN_TRA_par["SCpos_MAR_ini"]
    DYN_TRA_out["SCpos_SUN"] = DYN_TRA_par["SCpos_SUN_ini"]
    DYN_TRA_out["SCpos_SSB"] = DYN_TRA_par["SCpos_SSB_ini"]
    DYN_TRA_out["SCvel_SSB"] = DYN_TRA_par["SCvel_SSB_ini"]
    return dict(DYN_TRA_out)

# Runge-Kutta 4 integrator.
# TODO: Set this as a generic function in utils to be used on other modules
def rk4_step(dt, pos, vel, acc_func):
    k1_pos = vel
    k1_vel = acc_func(pos, vel)

    pos2 = pos + 0.5*dt*k1_pos
    vel2 = vel + 0.5*dt*k1_vel
    k2_pos = vel2
    k2_vel = acc_func(pos2, vel2)

    pos3 = pos + 0.5*dt*k2_pos
    vel3 = vel + 0.5*dt*k2_vel
    k3_pos = vel3
    k3_vel = acc_func(pos3, vel3)

    pos4 = pos + 0.5*dt*k3_pos
    vel4 = vel + 0.5*dt*k3_vel
    k4_pos = vel4
    k4_vel = acc_func(pos4, vel4)

    pos_next = pos + dt*(k1_pos + 2*k2_pos + 2*k3_pos * k4_pos)/6.0
    vel_next = vel + dt*(k1_vel + 2*k2_vel + 2*k3_vel * k4_vel)/6.0

    return pos_next, vel_next
