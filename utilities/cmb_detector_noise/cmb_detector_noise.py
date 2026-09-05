import numpy as np
import matplotlib.pyplot as plt

# Parameters
Tsys = 30    # [K] Overall system temperature with added noise
dt = 360     # [s] Integration time
f_min =  1e9 # [Hz] Minimum test bandwidth
f_max = 10e9 # [Hz] Maximum test bandwidth

# Compute minimum detectable temperature at the antenna
f = np.linspace(f_min, f_max, 10)
Ta = 2*Tsys/np.sqrt(f*dt) # [K] (rms)
Ta_uK = Ta*1e6 # [µK] (rms)

# For the detectable temperature range, compute the detectable velocity change over angle
theta = 0      # [rad] Angle between the velocity vector and the detector direction
T0_cmb = 2.73  # [K] CMB radiation mean temperature
c = 299792.458 # [km/s] Light speed
v_SSB = 369.8  # [km/s] SSB velocity
beta = v_SSB/c # [km/s] Normalized SSB velocity

T_factor = c * Ta/T0_cmb
theta = np.deg2rad(np.linspace(0, 89)) # [rad]
theta_factor = (1-beta*np.cos(theta))**2 * np.sqrt(1-beta**2) / (np.cos(theta) - beta)

plt.figure()
plt.plot(f, Ta_uK)
plt.xscale('log')
plt.xlabel("f [Hz]")
plt.ylabel("Ta [µK]")
plt.grid()

plt.figure()
plt.plot(np.rad2deg(theta), theta_factor)
plt.xlabel("theta [deg]")
plt.ylabel("theta factor")
plt.grid()

plt.figure()
plt.title(f"Velocity sensitivity for {dt} s of integration time")
# For each detectable temperature, plot the detectable velocity vs angle
for i, T in enumerate(T_factor):
    dv = T * theta_factor
    plt.plot(np.rad2deg(theta), dv, label=f"{f[i]/1e9:.1f} GHz")

plt.xlabel("Theta [deg]")
plt.ylabel("Velocity sensitivity [km/s]")
plt.legend()
plt.grid()

plt.show()
