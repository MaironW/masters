# Given:
# - A transfer date
# - Earth and Moon position relative to the Sun
# - Launch pad latitude
# - Spacecraft circular orbit elements
# - 3 reference stars
# Get:
# - The best time-fuel efficient transfer trajectory from Earth to the Moon
# - For each step of the orbit, compute the alpha angle between the Moon and the 3 reference stars

import numpy as np
import matplotlib.pyplot as plt
import astropy.units as u
from astropy.time import Time
from poliastro.bodies import Sun, Earth, Moon
from poliastro.twobody import Orbit
from poliastro.ephem import Ephem
from poliastro.plotting.static import StaticOrbitPlotter
from poliastro.maneuver import Maneuver


# Plot colors
c1 = (0.0,0.4,1.0)
c2 = (0.8,0.0,0.0)
c3 = (0.2,0.8,0.0)
c4 = 'm'

# Define direction of arbitrary Star relative to the Sun
# It must be inside the FOV of the Spacecraft to make sense
star_pos_sun_unit = [0.987, 0.101, 0.129] << u.one
star_pos_sun_unit = [1.000, 0.000, 0.000] << u.one
star_pos_sun_unit = [0.500, 0.500, 0.500] << u.one

# Launch data
launch_lat = 28.4 << u.deg # Latitude: Cape Canaveral
launch_az  = 90.0 << u.deg # Azimuth: east

# Time in Barycentric Dynamic Time
start_time = Time("2025-05-03", scale="tdb")

# Set Spacecraft Orbit parameters
h    = 400
a    = 6378 + h << u.km  # 400 km circular orbit
ecc  = 0.0001   << u.one
inc  = np.arccos(np.cos(launch_lat)*np.sin(launch_az))
raan = 0        << u.deg
argp = 0        << u.deg
tano = 0        << u.deg

# Generate Keplerian orbit at start time
keplerian_orbit = Orbit.from_classical(Earth, a, ecc, inc, raan, argp, tano, epoch=start_time)
kep_rv          = keplerian_orbit.rv()

# Get Moon orbit inclination
moon_ephem = Ephem.from_body(Moon, start_time, attractor=Earth)
moon_orbit = Orbit.from_ephem(Earth, moon_ephem, start_time)
moon_inc   = moon_orbit.inc.to(u.deg)

# Relative inclination between Satellite and Moon around Earth
del_inc = keplerian_orbit.inc.to(u.deg) - moon_inc

# Compute the transfer speed from the start orbit to the Moon inclined orbit
sat_v = np.linalg.norm(kep_rv[1])
dV_inc = sat_v * np.sqrt(2*(1-np.cos(del_inc)))

# From now on, we assume the satellite is orbiting on the same plane as the Moon
inc = moon_inc
keplerian_orbit = Orbit.from_classical(Earth, a, ecc, inc, raan, argp, tano, epoch=start_time)
kep_rv          = keplerian_orbit.rv()

# Create the transfer orbit
transfer_orbit_a   = (keplerian_orbit.r_p + moon_orbit.r_a)/2
transfer_orbit_ecc = 1 - keplerian_orbit.r_p / transfer_orbit_a
transfer_orbit = Orbit.from_classical(Earth,
                                      transfer_orbit_a,
                                      transfer_orbit_ecc,
                                      keplerian_orbit.inc,
                                      keplerian_orbit.raan,
                                      keplerian_orbit.argp,
                                      keplerian_orbit.nu,
                                      epoch=start_time)

# Plot true anomalies along the time (Propagate)
# This is poorly done, just for the sake of verification
propagation_time    = int(moon_orbit.period.to(u.hour).to_value())
transfer_orbit_tano = np.zeros(propagation_time)
moon_orbit_tano     = np.zeros(propagation_time)
time_range          = np.zeros(propagation_time)
for i in range(propagation_time):
    transfer_orbit_prop    = transfer_orbit.propagate(i << u.hour)
    moon_orbit_prop        = moon_orbit.propagate(i << u.hour)
    transfer_orbit_tano[i] = transfer_orbit_prop.nu.to(u.deg).to_value()%360
    moon_orbit_tano[i]     = moon_orbit_prop.nu.to(u.deg).to_value()%360
    time_range[i] = i

tano_compare    = abs(transfer_orbit_tano - moon_orbit_tano)
rendezvous_idx  = np.argmin(tano_compare)
rendezvous_time = time_range[rendezvous_idx] << u.hour
print(moon_orbit_tano[rendezvous_idx])
print(transfer_orbit_tano[rendezvous_idx])

transfer_orbit_rendezvous = transfer_orbit.propagate(rendezvous_time)
moon_orbit_rendezvous     = moon_orbit.propagate(rendezvous_time)

print(rendezvous_time)

plt.figure()
plt.plot(time_range, abs(transfer_orbit_tano - moon_orbit_tano))
plt.grid()

plt.figure()
plt.plot(time_range, transfer_orbit_tano,label='Spacecraft True Anomaly')
plt.plot(time_range, moon_orbit_tano,label='Moon True Anomaly')
plt.xlabel("t (days)")
plt.ylabel("f (deg)")
plt.grid()
plt.legend()

# Plot orbits
fig, ax = plt.subplots()
op = StaticOrbitPlotter(ax)
op.plot(keplerian_orbit, label='Spacecraft')
op.plot(moon_orbit, label='Moon')
op.plot(moon_orbit_rendezvous, label='Moon')
op.plot(transfer_orbit, label='Transfer (Homman)')
op.plot(transfer_orbit_rendezvous, label='Transfer (Homman)')
plt.legend()
plt.show()