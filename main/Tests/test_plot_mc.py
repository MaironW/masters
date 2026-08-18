# Sample Monte Carlo plot script
# Run from repository root with:
#   python3 -m Tests.test_plot_mc


from PPC   import PPC, PPC_MC

fig, ax, x, y_runs, y_mean = PPC_MC.plot_mc(log_path="Logs/test_mc", 
                                            x_var="DYN.DYN_TIME.time_SIM", y_var="NAV.NAV_EKF.x_est",
                                            xlabel="Time [s]", ylabel="y(t)",
                                            x_axis=None, y_axis=0,
                                            title="System Output y(t)")

PPC.show_plot()