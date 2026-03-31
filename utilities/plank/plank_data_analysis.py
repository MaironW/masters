import healpy as hp
import numpy as np
import matplotlib.pyplot as plt
import os

# File path
filename = "COM_CMB_IQU-smica_2048_R3.00_full.fits"

# Load Planck anisotropy map
cmb_map = hp.read_map(filename, field=0)

# OPTIONAL: downgrade resolution for speed
cmb_map = hp.ud_grade(cmb_map, nside_out=64)

# Constants
T0    = 2.725       # [K]
c     = 299792458.0 # [m/s]
v_mag = 371e3       # [m/s]

# Galactic velocity direction
l = np.deg2rad(263.85) # [rad]
b = np.deg2rad(48.25)  # [rad]

# Galactic velocity direction
v_hat = np.array([
    np.cos(b) * np.cos(l),
    np.cos(b) * np.sin(l),
    np.sin(b)
])
v = v_mag * v_hat

# Spacecraft velocity (approx)
v_sc_mag = 27e3 # [m/s]

# Spacecraft velocity Direction (Galactic)
l_sc = np.deg2rad(90.0)
b_sc = np.deg2rad(0.0)

v_sc_hat = np.array([
    np.cos(b_sc) * np.cos(l_sc),
    np.cos(b_sc) * np.sin(l_sc),
    np.sin(b_sc)
])

v_sc = v_sc_mag * v_sc_hat

v_total = v + v_sc

beta  = v_total / c
beta2 = np.dot(beta, beta)

# HEALPix geometry
nside = hp.get_nside(cmb_map)
npix  = hp.nside2npix(nside)

theta, phi = hp.pix2ang(nside, np.arange(npix))

n_hat = np.vstack([
    np.sin(theta) * np.cos(phi),
    np.sin(theta) * np.sin(phi),
    np.cos(theta)
]).T

# 1) ANISOTROPIES ONLY
anisotropy_map = cmb_map # Already delta T

# 2) DIPOLE ONLY (relativistic)
gamma      = np.sqrt(1 - beta2)
beta_dot_n = n_hat @ beta
dipole_map = T0 * gamma / (1 - beta_dot_n)

# Remove monopole to visualize dipole clearly
dipole_map = dipole_map - np.mean(dipole_map)

# 3) FULL MAP (anisotropy + dipole)
full_map = (T0 + cmb_map) * gamma / (1 - beta_dot_n)

# Diagnostics
print("Anisotropy RMS:", np.std(anisotropy_map))
print("Dipole peak-to-peak:", np.max(dipole_map) - np.min(dipole_map))
print("Full map mean:", np.mean(full_map))

# PLOTS

# Anisotropies
hp.mollview(
    anisotropy_map,
    title="Planck CMB Anisotropies (µK scale)",
    unit="K",
    cmap="coolwarm"
)
hp.graticule()

# Dipole
hp.mollview(
    dipole_map,
    title="CMB Dipole (Relativistic, velocity only)",
    unit="K",
    cmap="coolwarm"
)
hp.graticule()

# Full map
hp.mollview(
    full_map,
    title="CMB Full Sky (Dipole + Anisotropies)",
    unit="K",
    cmap="coolwarm"
)

hp.graticule()

##########################
# TIME SERIES SIMULATION #
##########################

# Spin parameters
spin_rate = 2 * np.pi / 60.0   # [rad/s] → 1 rotation per 60 s
t = np.linspace(0, 300, 2000)  # 5 minutes

# Sensor boresight in spacecraft frame (e.g. X-axis)
bore_sc = np.array([1.0, 0.0, 0.0])

T_time = []

for ti in t:
    # Rotation angle
    theta_spin = spin_rate * ti

    # Tilt spacecraft axis
    tilt_angle = np.deg2rad(30)
    Ry = np.array([
        [ np.cos(tilt_angle), 0, np.sin(tilt_angle)],
        [ 0,                 1, 0],
        [-np.sin(tilt_angle), 0, np.cos(tilt_angle)]
    ])

    # Rotation matrix around Z-axis
    Rz = np.array([
        [ np.cos(theta_spin), -np.sin(theta_spin), 0],
        [ np.sin(theta_spin),  np.cos(theta_spin), 0],
        [ 0,                  0,                 1]
    ])

    # Rotate boresight into inertial (Galactic) frame
    n_t = Ry @ (Rz @ bore_sc)

    # Convert to healpy angles
    theta_t = np.arccos(n_t[2])           # colatitude
    phi_t   = np.arctan2(n_t[1], n_t[0])  # longitude

    if phi_t < 0:
        phi_t += 2*np.pi

    # Sample map
    T_sample = hp.get_interp_val(full_map, theta_t, phi_t)

    T_time.append(T_sample)

T_time = np.array(T_time)
T_detrended = T_time - np.mean(T_time)

fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True)
plt.title("CMB Temperature vs Time (Spinning Spacecraft)")
ax1.plot(t, T_time)
ax2.plot(t, T_detrended)
ax2.set_xlabel("Time [s]")
ax1.set_ylabel("Temperature [K]")
ax2.set_ylabel("Temperature [K]")
ax1.grid()
ax2.grid()
plt.show()