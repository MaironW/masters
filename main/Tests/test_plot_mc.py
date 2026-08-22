# Sample Monte Carlo plot script
# Run from repository root with:
#   python3 -m Tests.test_plot_mc

import numpy as np

from PPC   import PPC, PPC_MC

runs, files = PPC_MC.load_mc_logs("Logs/test_mc")

def EKF_error(timeline):
    SCpos_SSB = timeline["DYN"]["DYN_TRA"]["SCpos_SSB"]
    x_est = timeline["NAV"]["NAV_EKF"]["x_est"]
    SCpos_SSB_est = x_est[:, 0:3]
    SCpos_SSB_err = SCpos_SSB_est - SCpos_SSB
    return SCpos_SSB_err

def covariance(timeline):
    P = timeline["NAV"]["NAV_EKF"]["P"]
    P_diag = P[:, 0:3, 0:3].diagonal(axis1=1, axis2=2)
    return np.sqrt(P_diag)

time_SIM, pos_err, pos_err_mean, pos_err_std = PPC_MC.mc_data(runs, postprocess=EKF_error, x_var="DYN.DYN_TIME.time_SIM",y_axis=0)
time_SIM_days = time_SIM/(24*3600)

time_SIM, cov, cov_mean, cov_std = PPC_MC.mc_data(runs, postprocess=covariance, x_var="DYN.DYN_TIME.time_SIM",y_axis=0)

fig, ax = PPC_MC.plot_mc(time_SIM_days, pos_err, pos_err_mean,
                         xlabel="Time [days]", ylabel="y(t)",
                         title="EKF Position Error [km]")

fig, ax = PPC.plot(time_SIM_days, pos_err_mean, label="Mean", fig=fig,ax=ax)
fig, ax = PPC.plot(time_SIM_days, pos_err_mean+3*pos_err_std, color=PPC.colors["red"], label="3σ", fig=fig,ax=ax)
fig, ax = PPC.plot(time_SIM_days, pos_err_mean-3*pos_err_std, color=PPC.colors["red"], fig=fig,ax=ax)

# The covariance does not change for each Monte Carlo run with the same covariance matrices
fig, ax = PPC.plot(time_SIM_days, +3*cov[0], color=PPC.colors["green"], label="3σ P", fig=fig,ax=ax)
fig, ax = PPC.plot(time_SIM_days, -3*cov[0], color=PPC.colors["green"], fig=fig,ax=ax)

# Note: filter overconfident, resulting in covariance way smaller than monte carlo std

PPC.show_plot()