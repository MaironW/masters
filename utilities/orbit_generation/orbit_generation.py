import numpy as np
from numpy import sin, cos, pi
import matplotlib.pyplot as plt
from astropy import units as u
from poliastro.bodies import Earth
from poliastro.twobody import Orbit
from poliastro.twobody.propagation import CowellPropagator
from poliastro.core.propagation import func_twobody
from poliastro.core.perturbations import J2_perturbation

#############
# FUNCTIONS #
#############

# Define perturbation model
def perturbation(t0, u_, k):
    du_kep = func_twobody(t0, u_, k)
    ax, ay, az = J2_perturbation(t0, u_, k, J2=Earth.J2.value, R=Earth.R.to(u.km).value)
    du_ad = np.array([0, 0, 0, ax, ay, az])
    return du_kep + du_ad

#################
# CONFIGURATION #
#################

# Set propagation time span
dt        = 1 << u.s
time_prop = 6 << u.h
time_span = np.arange(0, time_prop.to(u.s).value , dt.value) << u.s
n_steps   = len(time_span)

# Set Orbit parameters
a    = 6731     << u.km
ecc  = 0.0001   << u.one
inc  = 51       << u.deg
raan = 13       << u.deg
argp = 42       << u.deg
nu   = 23.33    << u.deg # true anomaly

##################
# GENERATE ORBIT #
##################

# Generate Keplerian Orbit
keplerian_orbit = Orbit.from_classical(Earth, a, ecc, inc, raan, argp, nu)

# Initialize state arrays
kep_rv   = np.zeros((n_steps,6))
per_rv   = np.zeros((n_steps,6))
kep_elem = np.zeros((n_steps,6))
per_elem = np.zeros((n_steps,6))

# Propagate non-perturbated (kep) and perturbated (per) orbits
for t_prop in range(n_steps):
    print(time_span[t_prop])
    kep_ephem = keplerian_orbit.propagate(time_span[t_prop])
    per_ephem = keplerian_orbit.propagate(time_span[t_prop], method=CowellPropagator(f=perturbation))

    kep_rv[t_prop,0:3] = kep_ephem.r.to(u.km)
    kep_rv[t_prop,3:6] = kep_ephem.v.to(u.km/u.s)

    per_rv[t_prop,0:3] = per_ephem.r.to(u.km)
    per_rv[t_prop,3:6] = per_ephem.v.to(u.km/u.s)

    kep_elem[t_prop, 0] = kep_ephem.a.to(u.km).value
    kep_elem[t_prop, 1] = kep_ephem.ecc.value
    kep_elem[t_prop, 2] = kep_ephem.inc.to(u.deg).value
    kep_elem[t_prop, 3] = kep_ephem.raan.to(u.deg).value
    kep_elem[t_prop, 4] = kep_ephem.argp.to(u.deg).value
    kep_elem[t_prop, 5] = kep_ephem.nu.to(u.deg).value

    per_elem[t_prop, 0] = per_ephem.a.to(u.km).value
    per_elem[t_prop, 1] = per_ephem.ecc.value
    per_elem[t_prop, 2] = per_ephem.inc.to(u.deg).value
    per_elem[t_prop, 3] = per_ephem.raan.to(u.deg).value
    per_elem[t_prop, 4] = per_ephem.argp.to(u.deg).value
    per_elem[t_prop, 5] = per_ephem.nu.to(u.deg).value

###############
# SAVE STATES #
###############
np.savetxt("kep_orbit.csv", kep_rv, delimiter=',')
np.savetxt("per_orbit.csv", per_rv, delimiter=',')

#########
# PLOTS #
#########

# Position
plt.figure('Position ECI (km)')
axs    = [0]*6
label_list = ["x (km)","y (km)","z (km)"]
for i in range(3):
    axs[i] = plt.subplot(3,1,i+1)
    axs[i].grid()
    axs[i].plot(time_span.to(u.h), kep_rv[:,i], label = 'Keplerian')
    axs[i].plot(time_span.to(u.h), per_rv[:,i], label = 'Perturbed')
    axs[i].legend()
    axs[i].set_ylabel(label_list[i])
axs[2].set_xlabel("Time (h)")

# Velocity
plt.figure('Velocity ECI (km/s)')
axs    = [0]*6
label_list = ["vx (km/s)","vy (km/s)","vz (km/s)"]
for i in range(3):
    axs[i] = plt.subplot(3,1,i+1)
    axs[i].grid()
    axs[i].plot(time_span.to(u.h), kep_rv[:,i+3], label = 'Keplerian')
    axs[i].plot(time_span.to(u.h), per_rv[:,i+3], label = 'Perturbed')
    axs[i].legend()
    axs[i].set_ylabel(label_list[i])
axs[2].set_xlabel("Time (h)")

# Orbital Elements
axs    = [0]*6
label_list = ["sma (km)","ecc","inc (deg)","raan (deg)","argp (deg)","tano (deg)"]
plt.figure("Orbital Elements")
for i in range(6):
    axs[i] = plt.subplot(6,1,i+1)
    axs[i].grid()
    axs[i].plot(time_span.to(u.h), kep_elem[:,i], label = 'Keplerian')
    axs[i].plot(time_span.to(u.h), per_elem[:,i], label = 'Perturbed')
    axs[i].legend()
    axs[i].set_ylabel(label_list[i])
axs[5].set_xlabel("Time (h)")

# 3D plot
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot(kep_rv[:, 0], kep_rv[:, 1], kep_rv[:, 2])
ax.plot(per_rv[:, 0], per_rv[:, 1], per_rv[:, 2])
ax.set_xlabel('X (km)')
ax.set_ylabel('Y (km)')
ax.set_zlabel('Z (km)')

# Make surface data
u_sur = np.linspace(0, 2*pi, 20)
v_sur = np.linspace(0, pi, 20)
x_sur = Earth.R.to(u.km).value * np.outer(cos(u_sur), sin(v_sur))
y_sur = Earth.R.to(u.km).value * np.outer(sin(u_sur), sin(v_sur))
z_sur = Earth.R.to(u.km).value * np.outer(np.ones(np.size(u_sur)), cos(v_sur))

# Plot the surface
ax.plot_surface(x_sur, y_sur, z_sur, color='C0', alpha=0.1)
ax.set_aspect('equal')

plt.show()
