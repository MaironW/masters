# Sample Monte Carlo plot script
# Run from repository root with:
#   python3 -m Tests.test_plot_mc

import numpy as np

from PPC   import PPC, PPC_MC

runs, files = PPC_MC.load_mc_logs("Logs/test_mc")

def EKF_pos_error(timeline):
    SCpos_SSB = timeline["DYN"]["DYN_TRA"]["SCpos_SSB"]
    x_est = timeline["NAV"]["NAV_EKF"]["x_est"]
    SCpos_SSB_est = x_est[:, 0:3]
    SCpos_SSB_err = SCpos_SSB_est - SCpos_SSB
    return SCpos_SSB_err

def EKF_vel_error(timeline):
    SCvel_SSB = timeline["DYN"]["DYN_TRA"]["SCvel_SSB"]
    x_est = timeline["NAV"]["NAV_EKF"]["x_est"]
    SCvel_SSB_est = x_est[:, 3:6]
    SCvel_SSB_err = SCvel_SSB_est - SCvel_SSB
    return SCvel_SSB_err

def pos_covariance(timeline):
    P = timeline["NAV"]["NAV_EKF"]["P"]
    P_diag = P[:, 0:3, 0:3].diagonal(axis1=1, axis2=2)
    return np.sqrt(P_diag)

def vel_covariance(timeline):
    P = timeline["NAV"]["NAV_EKF"]["P"]
    P_diag = P[:, 3:6, 3:6].diagonal(axis1=1, axis2=2)
    return np.sqrt(P_diag)

time_SIM, pos_err, pos_err_mean, pos_err_std = PPC_MC.mc_data(runs, postprocess=EKF_pos_error, x_var="DYN.DYN_TIME.time_SIM",y_axis=0)
time_SIM, vel_err, vel_err_mean, vel_err_std = PPC_MC.mc_data(runs, postprocess=EKF_vel_error, x_var="DYN.DYN_TIME.time_SIM",y_axis=0)
time_SIM_days = time_SIM/(24*3600)

time_SIM, pos_cov, pos_cov_mean, pos_cov_std = PPC_MC.mc_data(runs, postprocess=pos_covariance, x_var="DYN.DYN_TIME.time_SIM",y_axis=0)
time_SIM, vel_cov, vel_cov_mean, vel_cov_std = PPC_MC.mc_data(runs, postprocess=vel_covariance, x_var="DYN.DYN_TIME.time_SIM",y_axis=0)

fig, ax1 = PPC_MC.plot_mc(time_SIM_days, pos_err, pos_err_mean,
                         xlabel="Time [days]", ylabel="y(t)",
                         title="EKF Position Error [km]", subplot=(2,1,1))

fig, ax2 = PPC_MC.plot_mc(time_SIM_days, vel_err, vel_err_mean,
                         xlabel="Time [days]", ylabel="y(t)",
                         title="EKF Velocity Error [km/s]", fig=fig, subplot=(2,1,2))

# Position
fig, ax1 = PPC.plot(time_SIM_days, pos_err_mean, label="Mean", fig=fig,ax=ax1)
fig, ax1 = PPC.plot(time_SIM_days, pos_err_mean+3*pos_err_std, color=PPC.colors["red"], label="3σ", fig=fig,ax=ax1)
fig, ax1 = PPC.plot(time_SIM_days, pos_err_mean-3*pos_err_std, color=PPC.colors["red"], fig=fig,ax=ax1)
# The covariance does not change for each Monte Carlo run with the same covariance matrices
fig, ax1 = PPC.plot(time_SIM_days, +3*pos_cov[0], color=PPC.colors["green"], label="3σ P", fig=fig,ax=ax1)
fig, ax1 = PPC.plot(time_SIM_days, -3*pos_cov[0], color=PPC.colors["green"], fig=fig,ax=ax1)

# Velocity
fig, ax2 = PPC.plot(time_SIM_days, vel_err_mean, label="Mean", fig=fig,ax=ax2)
fig, ax2 = PPC.plot(time_SIM_days, vel_err_mean+3*vel_err_std, color=PPC.colors["red"], label="3σ", fig=fig,ax=ax2)
fig, ax2 = PPC.plot(time_SIM_days, vel_err_mean-3*vel_err_std, color=PPC.colors["red"], fig=fig,ax=ax2)
# The covariance does not change for each Monte Carlo run with the same covariance matrices
fig, ax2 = PPC.plot(time_SIM_days, +3*vel_cov[0], color=PPC.colors["green"], label="3σ P", fig=fig,ax=ax2)
fig, ax2 = PPC.plot(time_SIM_days, -3*vel_cov[0], color=PPC.colors["green"], fig=fig,ax=ax2)

# Note: filter overconfident, resulting in covariance way smaller than monte carlo std

PPC.show_plot()