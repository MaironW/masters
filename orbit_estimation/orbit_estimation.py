import numpy as np
from numpy import pi
import matplotlib.pyplot as plt
from scipy.stats import multivariate_normal
from kalman_filter import kalman_filter
from particle_filter import particle_filter
from astropy import units as units
from orbit_conversion import elementsfromstate

# Plot colors
c1 = (0.0,0.4,1.0)
c2 = (0.8,0.0,0.0)
c3 = (0.2,0.8,0.0)
c4 = 'm'

####################
# DEFINE CONSTANTS #
####################
mu = 3.986e5 # [km^3/s^2] Earth gravitational parameter

#####################
# LOAD MEASUREMENTS #
#####################

kep_rv = np.loadtxt("../orbit_generation/kep_orbit.csv",delimiter=',')[:10000]
per_rv = np.loadtxt("../orbit_generation/per_orbit.csv",delimiter=',')[:10000]

# Select between the Keplerian or Pertubed orbit
rv = per_rv

###################
# SET FILTER TIME #
###################

dt     = 1 # [s] data time step
mes_dt = 60 # [s] measurment time step
n_iters, n_states = rv.shape # [samples]
time_span = np.arange(0, n_iters*dt, dt)

##############
# STATISTICS #
##############

# Measurement Covariance (R)
pos_var = 1.0   # ± [km]
vel_var = 10e-3 # ± [km/s]

R = np.array([
    [pos_var**2, 0, 0, 0, 0, 0],
    [0, vel_var**2, 0, 0, 0, 0],
    [0, 0, pos_var**2, 0, 0, 0],
    [0, 0, 0, vel_var**2, 0, 0],
    [0, 0, 0, 0, pos_var**2, 0],
    [0, 0, 0, 0, 0, vel_var**2],
])

# Process Noise Covariance (Q)
# TODO: Just placeholder, in the order of J2
acc_var = 1e-3 # ± [km/s^2]
Q = acc_var*acc_var*np.array([[0.25*dt**4, 0.5*dt**3,          0,         0,          0,         0],
                              [0.5*dt**3,      dt**2,          0,         0,          0,         0],
                              [0,                  0, 0.25*dt**4, 0.5*dt**3,          0,         0],
                              [0,                  0,  0.5*dt**3,     dt**2,          0,         0],
                              [0,                  0,          0,         0, 0.25*dt**4, 0.5*dt**3],
                              [0,                  0,          0,         0,  0.5*dt**3,     dt**2]])

# Create distribuitions for the models
# allow_singular=True must be used to allow ill_conditioned matrix
# TODO: learn how to make the matrix Q not ill conditioned
measurement_distribuition = multivariate_normal(mean=[0]*n_states, cov=R)
process_distribuition     = multivariate_normal(mean=[0]*n_states, cov=Q, allow_singular=True)

measurement_noise = measurement_distribuition.rvs(size=n_iters)

#####################
# STATE-SPACE MODEL #
#####################

r = rv[:,0:3]
v = rv[:,3:6]

# Simulate measurement delay
waiting_time = int(mes_dt/dt)
rx_new = np.zeros(n_iters)
ry_new = np.zeros(n_iters)
rz_new = np.zeros(n_iters)
vx_new = np.zeros(n_iters)
vy_new = np.zeros(n_iters)
vz_new = np.zeros(n_iters)
for i in range(0, n_iters, waiting_time):
    rx_new[i:i+waiting_time] = r[i,0] #+ measurement_noise[i,0]
    ry_new[i:i+waiting_time] = r[i,1] #+ measurement_noise[i,2]
    rz_new[i:i+waiting_time] = r[i,2] #+ measurement_noise[i,4]
    vx_new[i:i+waiting_time] = v[i,0] #+ measurement_noise[i,1]
    vy_new[i:i+waiting_time] = v[i,1] #+ measurement_noise[i,3]
    vz_new[i:i+waiting_time] = v[i,2] #+ measurement_noise[i,5]
rx_mes = np.copy(rx_new)
ry_mes = np.copy(ry_new)
rz_mes = np.copy(rz_new)
vx_mes = np.copy(vx_new)
vy_mes = np.copy(vy_new)
vz_mes = np.copy(vz_new)

# Linear Dynamic matrix (just for Kalman Filter)
kf_A = np.array([
    [1, dt, 0,  0, 0,  0],
    [0,  1, 0,  0, 0,  0],
    [0,  0, 1, dt, 0,  0],
    [0,  0, 0,  1, 0,  0],
    [0,  0, 0,  0, 1, dt],
    [0,  0, 0,  0, 0,  1]
    ])

# Dynamic matrix for Particle Filter
# It is actually a function of the position state
# Because of this, the linear Kalman Filter is not applicable
def F(states, dt):
    mu = 3.986e5 # [km^3/s^2] Earth gravitational parameter
    n_states, n_particles = states.shape
    for j in range(n_particles):
        rx  = states[0,j]
        ry  = states[2,j]
        rz  = states[4,j]
        vx  = states[1,j]
        vy  = states[3,j]
        vz  = states[5,j]
        R_3 = (rx*rx + ry*ry + rz*rz)**1.5
        tmp = -mu/R_3
        ax  = tmp*rx
        ay  = tmp*ry
        az  = tmp*rz

        # propagate forward
        vx_new = vx + ax * dt
        vy_new = vy + ay * dt
        vz_new = vz + az * dt
        rx_new = rx + vx_new * dt
        ry_new = ry + vy_new * dt
        rz_new = rz + vz_new * dt

        states[:,j] = np.array([rx_new, vx_new, ry_new, vy_new, rz_new, vz_new])
    return states

# Control matrix (unused)
B = np.zeros((n_states, 3))

# Measurement matrix
H = np.eye(n_states)

# Measurement vector
# Rearange measurements so the state is actually [rx, vx, ry, vy, rz, vz]
# Z = np.vstack((r[:,0],v[:,0],r[:,1],v[:,1],r[:,2],v[:,2])).T
Z = np.vstack((rx_mes,vx_mes,ry_mes,vy_mes,rz_mes,vz_mes)).T

# Control input vector (unused)
u = np.zeros((3,n_iters)).T

###############
# RUN FILTERS #
###############

KF_est = kalman_filter(kf_A,B,H,Z,u,R,Q)
PF_est, PF_est_list, PF_weight_list = particle_filter(F,B,Z,u,R,Q,dt,n_particles=500)

KF_rv = np.vstack((KF_est[:,0], KF_est[:,2], KF_est[:,4], KF_est[:,1], KF_est[:,3], KF_est[:,5])).T
PF_rv = np.vstack((PF_est[:,0], PF_est[:,2], PF_est[:,4], PF_est[:,1], PF_est[:,3], PF_est[:,5])).T

################
# GET RESIDUES #
################

kep_r = kep_rv[:,0:3]
kep_v = kep_rv[:,3:6]

# Kalman Filter Errors
dx_est_err = abs(KF_rv[:,0]-r[:,0])
dy_est_err = abs(KF_rv[:,1]-r[:,1])
dz_est_err = abs(KF_rv[:,2]-r[:,2])
KF_pos_err = (dx_est_err**2 + dy_est_err**2 + dz_est_err**2)**0.5
vx_est_err = abs(KF_rv[:,3]-v[:,0])
vy_est_err = abs(KF_rv[:,4]-v[:,1])
vz_est_err = abs(KF_rv[:,5]-v[:,2])
KF_vel_err = (vx_est_err**2 + vy_est_err**2 + vz_est_err**2)**0.5

dx_est_err = abs(KF_rv[:,0]-kep_r[:,0])
dy_est_err = abs(KF_rv[:,1]-kep_r[:,1])
dz_est_err = abs(KF_rv[:,2]-kep_r[:,2])
KF_kep_pos_err = (dx_est_err**2 + dy_est_err**2 + dz_est_err**2)**0.5
vx_est_err = abs(KF_rv[:,3]-kep_v[:,0])
vy_est_err = abs(KF_rv[:,4]-kep_v[:,1])
vz_est_err = abs(KF_rv[:,5]-kep_v[:,2])
KF_kep_vel_err = (vx_est_err**2 + vy_est_err**2 + vz_est_err**2)**0.5

# Particle Filter Errors
dx_est_err = abs(PF_rv[:,0]-r[:,0])
dy_est_err = abs(PF_rv[:,1]-r[:,1])
dz_est_err = abs(PF_rv[:,2]-r[:,2])
PF_pos_err = (dx_est_err**2 + dy_est_err**2 + dz_est_err**2)**0.5
vx_est_err = abs(PF_rv[:,3]-v[:,0])
vy_est_err = abs(PF_rv[:,4]-v[:,1])
vz_est_err = abs(PF_rv[:,5]-v[:,2])
PF_vel_err = (vx_est_err**2 + vy_est_err**2 + vz_est_err**2)**0.5

dx_est_err = abs(PF_rv[:,0]-kep_r[:,0])
dy_est_err = abs(PF_rv[:,1]-kep_r[:,1])
dz_est_err = abs(PF_rv[:,2]-kep_r[:,2])
PF_kep_pos_err = (dx_est_err**2 + dy_est_err**2 + dz_est_err**2)**0.5
vx_est_err = abs(PF_rv[:,3]-kep_v[:,0])
vy_est_err = abs(PF_rv[:,4]-kep_v[:,1])
vz_est_err = abs(PF_rv[:,5]-kep_v[:,2])
PF_kep_vel_err = (vx_est_err**2 + vy_est_err**2 + vz_est_err**2)**0.5

##########################
# COMPUTE ORBIT ELEMENTS #
##########################
kep_r = kep_rv[:,0:3] << units.km
kep_v = kep_rv[:,3:6] << units.km/units.s
per_r = per_rv[:,0:3] << units.km
per_v = per_rv[:,3:6] << units.km/units.s
KF_r = KF_rv[:,0:3] << units.km
KF_v = KF_rv[:,3:6] << units.km/units.s
PF_r = PF_rv[:,0:3] << units.km
PF_v = PF_rv[:,3:6] << units.km/units.s

kep_elem = np.zeros((n_iters,6))
per_elem = np.zeros((n_iters,6))
KF_elem = np.zeros((n_iters,6))
PF_elem = np.zeros((n_iters,6))

kep_a, kep_ecc, kep_inc, kep_argp, kep_raan, kep_nu = elementsfromstate(kep_rv[:,0:3],kep_rv[:,3:6],mu)
per_a, per_ecc, per_inc, per_argp, per_raan, per_nu = elementsfromstate(per_rv[:,0:3],per_rv[:,3:6],mu)
KF_a, KF_ecc, KF_inc, KF_argp, KF_raan, KF_nu = elementsfromstate(KF_rv[:,0:3],KF_rv[:,3:6],mu)
PF_a, PF_ecc, PF_inc, PF_argp, PF_raan, PF_nu = elementsfromstate(PF_rv[:,0:3],PF_rv[:,3:6],mu)

kep_elem[:,0] = (kep_a << units.km).value
kep_elem[:,1] = kep_ecc
kep_elem[:,2] = (kep_inc << units.rad).to(units.deg).value
kep_elem[:,3] = (kep_raan << units.rad).to(units.deg).value
kep_elem[:,4] = (kep_argp << units.rad).to(units.deg).value
kep_elem[:,5] = (kep_nu << units.rad).to(units.deg).value
kep_elem[:,4] = ((kep_elem[:,4]+180)%360)-180
kep_elem[:,5] = ((kep_elem[:,5]+180)%360)-180

per_elem[:,0] = (per_a << units.km).value
per_elem[:,1] = per_ecc
per_elem[:,2] = (per_inc << units.rad).to(units.deg).value
per_elem[:,3] = (per_raan << units.rad).to(units.deg).value
per_elem[:,4] = (per_argp << units.rad).to(units.deg).value
per_elem[:,5] = (per_nu << units.rad).to(units.deg).value
per_elem[:,4] = ((per_elem[:,4]+180)%360)-180
per_elem[:,5] = ((per_elem[:,5]+180)%360)-180

KF_elem[:,0] = (KF_a << units.km).value
KF_elem[:,1] = KF_ecc
KF_elem[:,2] = (KF_inc << units.rad).to(units.deg).value
KF_elem[:,3] = (KF_raan << units.rad).to(units.deg).value
KF_elem[:,4] = (KF_argp << units.rad).to(units.deg).value
KF_elem[:,5] = (KF_nu << units.rad).to(units.deg).value
KF_elem[:,4] = ((KF_elem[:,4]+180)%360)-180
KF_elem[:,5] = ((KF_elem[:,5]+180)%360)-180

PF_elem[:,0] = (PF_a << units.km).value
PF_elem[:,1] = PF_ecc
PF_elem[:,2] = (PF_inc << units.rad).to(units.deg).value
PF_elem[:,3] = (PF_raan << units.rad).to(units.deg).value
PF_elem[:,4] = (PF_argp << units.rad).to(units.deg).value
PF_elem[:,5] = (PF_nu << units.rad).to(units.deg).value
PF_elem[:,4] = ((PF_elem[:,4]+180)%360)-180
PF_elem[:,5] = ((PF_elem[:,5]+180)%360)-180

################
# PLOT RESULTS #
################

# Position
plt.figure('Position ECI (km)')
axs    = [0]*6
label_list = ["x (km)","y (km)","z (km)"]
for i in range(3):
    axs[i] = plt.subplot(3,1,i+1)
    axs[i].grid()
    axs[i].plot(time_span/3600, kep_rv[:,i], color=c1, label = 'Keplerian')
    axs[i].plot(time_span/3600, per_rv[:,i], color=c2, label = 'Perturbed')
    axs[i].plot(time_span/3600, KF_rv[:,i], color=c3, label = 'Kalman Filter')
    axs[i].plot(time_span/3600, PF_rv[:,i], color=c4, label = 'Particle Filter')
    axs[i].legend()
    axs[i].set_ylabel(label_list[i])
axs[2].set_xlabel("Time (h)")

# Estimation Error
plt.figure("Estimation Error")
ax1 = plt.subplot(2,1,1)
ax1.grid()
ax1.plot(time_span/3600, KF_pos_err,  label='KF',  color=c3)
ax1.plot(time_span/3600, KF_kep_pos_err, '--',  label='KF Kep',  color=c3)
ax1.plot(time_span/3600, PF_pos_err,  label='PF',  color=c4)
ax1.plot(time_span/3600, PF_kep_pos_err, '--',  label='PF Kep',  color=c4)
ax1.set_ylabel("Pos (km)")
ax1.legend()

ax2 = plt.subplot(2,1,2)
ax2.grid()
ax2.plot(time_span/3600, KF_vel_err,  label='KF',  color=c3)
ax2.plot(time_span/3600, KF_kep_vel_err, '--',  label='KF Kep',  color=c3)
ax2.plot(time_span/3600, PF_vel_err,  label='PF',  color=c4)
ax2.plot(time_span/3600, PF_kep_vel_err, '--',  label='PF Kep',  color=c4)
ax2.set_ylabel("Vel (km/s)")
ax2.set_xlabel("Time (h)")
ax2.legend()

# Orbital Elements
axs = [0]*6
label_list = ["sma (km)","ecc","inc (deg)","raan (deg)","argp (deg)","tano (deg)"]
plt.figure("Orbital Elements")
for i in range(6):
    axs[i] = plt.subplot(6,1,i+1)
    axs[i].grid()
    axs[i].plot(time_span/3600, kep_elem[:,i], label = 'Keplerian', color=c1)
    axs[i].plot(time_span/3600, per_elem[:,i], label = 'Perturbed',color=c2)
    axs[i].plot(time_span/3600, KF_elem[:,i],  label = 'Kalman Filter',color=c3)
    axs[i].plot(time_span/3600, PF_elem[:,i],  label = 'Particle Filter',color=c4)
    axs[i].legend()
    axs[i].set_ylabel(label_list[i])
axs[5].set_xlabel("Time (h)")

# 3D plot
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot(kep_rv[:, 0], kep_rv[:, 1], kep_rv[:, 2], color=c1, label='Keplerian')
ax.plot(per_rv[:, 0], per_rv[:, 1], per_rv[:, 2], color=c2, label='Perturbed')
ax.plot(KF_rv[:, 0], KF_rv[:, 1], KF_rv[:, 2], color=c3, label='Kalman Filter')
ax.plot(PF_rv[:, 0], PF_rv[:, 1], PF_rv[:, 2], color=c4, label='Particle Filter')
ax.set_xlabel('X (km)')
ax.set_ylabel('Y (km)')
ax.set_zlabel('Z (km)')
ax.legend()

ax.set_aspect('equal')

plt.show()