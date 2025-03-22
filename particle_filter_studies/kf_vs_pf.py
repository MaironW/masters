# Compare a Kalman Filter and a Particle Filter for the same dynamic model
# Mairon de Souza Wolniewicz
# first version: 15/03/2025

import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import multivariate_normal
from numpy.random import randn

# Filter implementation functions
from kalman_filter   import kalman_filter
from particle_filter import particle_filter

# Plot colors
c1 = (0.0,0.4,1.0)
c2 = (0.8,0.0,0.0)
c3 = (0.2,0.8,0.0)

################################
# SET SIMULATION + FILTER TIME #
################################
dt      = 0.1 # [s] sim + filter time step
gps_dt  = 1.0 # [s] GPS measurement update step
n       = 8   #     number of cycles
t       = np.arange(0,n*2*np.pi,dt)
n_iters = len(t)

##############
# STATISTICS #
##############

# Measurement Covariance (R)
pos_var = 2.5 # ± [m]
vel_var = 0.1 # ± [m/s]
R = np.array([[pos_var*pos_var,0,0,0],
              [0,vel_var*vel_var,0,0],
              [0,0,pos_var*pos_var,0],
              [0,0,0,vel_var*vel_var]])

# Process Noise Covariance (Q)
acc_var = 0.32 # ± [m/s^2]
Q = acc_var*acc_var*np.array([[0.25*dt**4, 0.5*dt**3,          0,         0],
                              [0.5*dt**3,      dt**2,          0,         0],
                              [0,                  0, 0.25*dt**4, 0.5*dt**3],
                              [0,                  0,  0.5*dt**3,     dt**2]])

# Create distribuitions for the models
# allow_singular=True must be used to allow ill_conditioned matrix
# TODO: learn how to make the matrix Q not ill conditioned
measurement_distribuition = multivariate_normal(mean=[0,0,0,0], cov=R)
process_distribuition     = multivariate_normal(mean=[0,0,0,0], cov=Q, allow_singular=True)

##################
# DYNAMIC SYSTEM #
##################

# "Eight" trajectory
pathsize = 15 # [m] trajectory amplitude
# Position
dx = pathsize*np.cos(0.125*t) # [m]
dy = pathsize*np.sin(0.25*t)  # [m]
# Velocity
vx = -0.125*pathsize*np.sin(0.125*t) # [m/s]
vy = 0.25*pathsize*np.cos(0.25*t)    # [m/s]
# Acceleration
ax = -0.015625*pathsize*np.cos(0.125*t) # [m/s^2]
ay = -0.0625*pathsize*np.sin(0.25*t)    # [m/s^2]

measurement_noise = measurement_distribuition.rvs(size=n_iters)
process_noise     = process_distribuition.rvs(size=n_iters)

# Apply noise if desired
dx_noise = measurement_noise[:,0]
vx_noise = measurement_noise[:,1]
ax_noise = acc_var*randn(*ax.shape)
dy_noise = measurement_noise[:,2]
vy_noise = measurement_noise[:,3]
ay_noise = acc_var*randn(*ay.shape)

dx_gps = dx #+ dx_noise
dy_gps = dy #+ dy_noise
vx_gps = vx #+ vx_noise
vy_gps = vy #+ vy_noise
ax     = ax + ax_noise
ay     = ay + ay_noise

# Simulate GPS delay
waiting_time = int(gps_dt/dt)
dx_new = np.zeros(n_iters)
dy_new = np.zeros(n_iters)
vx_new = np.zeros(n_iters)
vy_new = np.zeros(n_iters)
for i in range(0, n_iters, waiting_time):
    dx_new[i:i+waiting_time]=dx_gps[i]
    dy_new[i:i+waiting_time]=dy_gps[i]
    vx_new[i:i+waiting_time]=vx_gps[i]
    vy_new[i:i+waiting_time]=vy_gps[i]
dx_gps = np.copy(dx_new)
dy_gps = np.copy(dy_new)
vx_gps = np.copy(vx_new)
vy_gps = np.copy(vy_new)

# Simulate invalid GPS data for a while
# START = int(len(t)*0.2)
# STOP = int(len(t)*0.28)
# dx_gps[START:STOP] = dx_gps[START]
# dy_gps[START:STOP] = dy_gps[START]
# vx_gps[START:STOP] = vx_gps[START]
# vy_gps[START:STOP] = vy_gps[START]

#####################
# STATE-SPACE MODEL #
#####################

# Dynamic matrix
A = np.array([
    [1, dt, 0,  0],
    [0,  1, 0,  0],
    [0,  0, 1, dt],
    [0,  0, 0,  1]
    ])

# Control matrix
B = np.array([
    [0.5*dt**2, 0],
    [dt,        0],
    [0, 0.5*dt**2],
    [0,        dt]
    ])

# Measurement matrix
H = np.array([
    [1, 0, 0, 0],
    [0, 1, 0, 0],
    [0, 0, 1, 0],
    [0, 0, 0, 1]
    ])

# Measurement vector
Z = np.array([dx_gps,vx_gps, dy_gps, vy_gps]).T

# Control input vector
u = np.array([ax,ay]).T

###############
# RUN FILTERS #
###############

KF_est = kalman_filter(A,B,H,Z,u,R,Q)
PF_est, PF_est_list, PF_weight_list = particle_filter(A,B,Z,u,R,Q,n_particles=300)

##################
# COMPUTE ERRORS #
##################

# Kalman Filter Errors
dx_est_err = abs(KF_est[:,0]-dx)
dy_est_err = abs(KF_est[:,2]-dy)
KF_pos_err = (dx_est_err**2 + dy_est_err**2)**0.5
vx_est_err = abs(KF_est[:,1]-vx)
vy_est_err = abs(KF_est[:,3]-vy)
KF_vel_err = (vx_est_err**2 + vy_est_err**2)**0.5

# Particle Filter Errors
dx_est_err = abs(PF_est[:,0]-dx)
dy_est_err = abs(PF_est[:,2]-dy)
PF_pos_err = (dx_est_err**2 + dy_est_err**2)**0.5
vx_est_err = abs(PF_est[:,1]-vx)
vy_est_err = abs(PF_est[:,3]-vy)
PF_vel_err = (vx_est_err**2 + vy_est_err**2)**0.5

################
# PLOT RESULTS #
################

# X-axis
plt.figure("X-axis")
ax1 = plt.subplot(3,1,1)
ax1.grid()
ax1.plot(t,          dx, label='true state',  color='k')
ax1.plot(t,      dx_gps, label='measurement', color=c1)
ax1.plot(t, KF_est[:,0], '--', label='KF',    color=c2)
ax1.plot(t, PF_est[:,0], '--', label='PF',    color=c3)
ax1.set_ylabel("Pos (m)")
ax1.legend()

ax2 = plt.subplot(3,1,2, sharex=ax1)
ax2.grid()
ax2.plot(t,          vx, label='true state',  color='k')
ax2.plot(t,      vx_gps, label='measurement', color=c1)
ax2.plot(t, KF_est[:,1], '--', label='KF',    color=c2)
ax2.plot(t, PF_est[:,1], '--', label='PF',    color=c3)
ax2.set_ylabel("Vel (m/s)")
ax2.legend()

ax3 = plt.subplot(3,1,3, sharex=ax1)
ax3.grid()
ax3.plot(t, ax, label='control input', color='k')
ax3.set_ylabel("Acc (m/s^2)")
ax3.set_xlabel("Time (s)")
ax3.legend()

# Y-axis
plt.figure("Y-axis")
ax1 = plt.subplot(3,1,1)
ax1.grid()
ax1.plot(t,          dy, label='true state',  color='k')
ax1.plot(t,      dy_gps, label='measurement', color=c1)
ax1.plot(t, KF_est[:,2], '--', label='KF',    color=c2)
ax1.plot(t, PF_est[:,2], '--', label='PF',    color=c3)
ax1.set_ylabel("Pos (m)")
ax1.legend()

ax2 = plt.subplot(3,1,2, sharex=ax1)
ax2.plot(t,          vy, label='true state',  color='k')
ax2.plot(t,      vy_gps, label='measurement', color=c1)
ax2.plot(t, KF_est[:,3], '--', label='KF',    color=c2)
ax2.plot(t, PF_est[:,3], '--', label='PF',    color=c3)
ax2.grid()
ax2.set_ylabel("Vel (m/s)")
ax2.legend()

ax3 = plt.subplot(3,1,3, sharex=ax1)
ax3.grid()
ax3.plot(t, ay, label='control input', color='k')
ax3.set_ylabel("Acc (m/s^2)")
ax3.set_xlabel("Time (s)")
ax3.legend()

# Estimation Error
plt.figure("Estimation Error")
ax1 = plt.subplot(2,1,1)
ax1.grid()
ax1.plot(t,KF_pos_err, label='KF', color=c2)
ax1.plot(t,PF_pos_err, label='PF', color=c3)
ax1.set_ylabel("Pos (m)")
ax1.legend()

ax2 = plt.subplot(2,1,2)
ax2.grid()
ax2.plot(t,KF_vel_err, label='KF', color=c2)
ax2.plot(t,PF_vel_err, label='PF', color=c3)
ax2.set_ylabel("Vel (m/s)")
ax2.set_xlabel("Time (s)")
ax2.legend()

# 2D Path
plt.figure("2D Path")
plt.grid()
plt.plot(dx_gps,dy_gps,'o-',color=c1,alpha=0.7,label='measurements')

plt.plot(KF_est[:,0],  KF_est[:,2],  '-', color=c2, label='KF')
plt.plot(KF_est[0,0],  KF_est[0,2],  'x', color=c2)
plt.plot(KF_est[-1,0], KF_est[-1,2], 'x', color=c2)

plt.plot(PF_est[:,0],  PF_est[:,2],  '-', color=c3, label='PF')
plt.plot(PF_est[0,0],  PF_est[0,2],  'x', color=c3)
plt.plot(PF_est[-1,0], PF_est[-1,2], 'x', color=c3)
plt.plot(dx,dy,color='k',label='true state')
plt.plot(dx[-1],dy[-1],'+k')

plt.xlabel("X (m)")
plt.ylabel("Y (m)")
plt.axis('equal')
plt.legend()
plt.plot(dx[0],dy[0],'o',color='orange')

plt.show()