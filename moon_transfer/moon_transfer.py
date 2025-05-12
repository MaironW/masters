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

from astropy.coordinates import solar_system_ephemeris, get_body_barycentric
from poliastro.bodies import Sun, Earth, Moon
from astropy.coordinates import CartesianRepresentation

# Plot colors
c1 = (0.0,0.4,1.0)
c2 = (0.8,0.0,0.0)
c3 = (0.2,0.8,0.0)
c4 = 'm'

# Define direction of arbitrary Star relative to the Sun
# It must be inside the FOV of the Spacecraft to make sense
star1_pos_sun_unit = [0.987, 0.101, 0.129] << u.one
star2_pos_sun_unit = [1.000, 0.000, 0.000] << u.one
star3_pos_sun_unit = [0.500, 0.500, 0.500] << u.one

# Launch data
launch_lat = 28.4 << u.deg # Latitude: Cape Canaveral
launch_az  = 90.0 << u.deg # Azimuth: east

# Time in Barycentric Dynamic Time
start_time = Time("2025-05-04T00:00:00", scale="tdb")

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

# Compute the transfer speed from the start orbit to the Moon inclined orbit
del_inc = keplerian_orbit.inc.to(u.rad) - moon_orbit.inc.to(u.rad)
sat_v = np.linalg.norm(kep_rv[1])
dV_inc = sat_v * np.sqrt(2*(1-np.cos(del_inc)))

# Compute the transfer speed to keep the orbit with the same argument of perigee of the Moon
del_argp = keplerian_orbit.argp.to(u.rad) - moon_orbit.argp.to(u.rad)
dV_argp = 2 * sat_v * np.sin(del_argp/2)

# From now on, we assume the satellite is orbiting on the same plane as the Moon
inc  = moon_orbit.inc.to(u.deg)
argp = moon_orbit.argp.to(u.deg)
raan = moon_orbit.raan.to(u.deg)

keplerian_orbit = Orbit.from_classical(Earth, a, ecc, inc, raan, argp, tano, epoch=start_time)

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
for i in range(propagation_time):# Compute the visible angle between the reference Star and the Moon
    transfer_orbit_prop    = transfer_orbit.propagate(i << u.hour)
    moon_orbit_prop        = moon_orbit.propagate(i << u.hour)
    transfer_orbit_tano[i] = transfer_orbit_prop.nu.to(u.deg).to_value()%360
    moon_orbit_tano[i]     = moon_orbit_prop.nu.to(u.deg).to_value()%360
    time_range[i] = i

tano_compare = abs(transfer_orbit_tano - moon_orbit_tano)
arrival_idx  = np.argmin(tano_compare)
arrival_time = time_range[arrival_idx] << u.day
arrival_tano = transfer_orbit_tano[arrival_idx] << u.deg

# This is not working very well. It would be best to do it algebraic
transfer_orbit_arrival = transfer_orbit.propagate(arrival_time)
moon_orbit_arrival     = moon_orbit.propagate(arrival_time)

# Plots
plt.figure()
plt.plot(time_range, abs(transfer_orbit_tano - moon_orbit_tano))
plt.xlabel("t (hours)")
plt.ylabel("f (deg)")
plt.grid()

plt.figure()
plt.plot(time_range, transfer_orbit_tano,label='Spacecraft True Anomaly')
plt.plot(time_range, moon_orbit_tano,label='Moon True Anomaly')
plt.xlabel("t (hours)")
plt.ylabel("f (deg)")
plt.grid()
plt.legend()

# Plot orbits
fig, ax = plt.subplots()
op = StaticOrbitPlotter(ax)
op.plot(keplerian_orbit, label='Spacecraft')
op.plot(moon_orbit, label='Moon')
op.plot(transfer_orbit, label='Transfer (Homman)')
op.plot(moon_orbit_arrival, label='Moon')
op.plot(transfer_orbit_arrival, label='Transfer (Homman)')
plt.legend()

###########
# CeleNav #
###########

# However, nothing will stop us againt imagining that there are Star Trackers pointing everywhere
# And that the satellite is aways looking to the Moon and the three reference stars

dt = 200
arrival_time = start_time + arrival_time
time_range = np.linspace(start_time,arrival_time,dt)
n_steps = len(time_range)
alpha1 = np.zeros(n_steps) << u.rad
alpha2 = np.zeros(n_steps) << u.rad
alpha3 = np.zeros(n_steps) << u.rad

with solar_system_ephemeris.set("builtin"):
    # Get Earth's position relative to the barycentric referential
    sun_pos_icrs   = get_body_barycentric("Sun", start_time)
    
    for i in range(n_steps):
        earth_pos_icrs = get_body_barycentric("Earth", time_range[i])
        moon_pos_icrs  = get_body_barycentric("Moon",  time_range[i])

        # Convert positions to the Sun reference frame
        earth_pos_sun = earth_pos_icrs - sun_pos_icrs
        moon_pos_sun = moon_pos_icrs   - sun_pos_icrs

        # Convert to cartesian representation
        earth_pos_sun = earth_pos_sun.represent_as(CartesianRepresentation).xyz.to(u.AU)
        moon_pos_sun  = moon_pos_sun.represent_as(CartesianRepresentation).xyz.to(u.AU)

        # Get position of Spacecraft relative to the Sun (iterate here)
        transfer_orbit = transfer_orbit.propagate(time_range[i])
        transfer_rv   = transfer_orbit.rv()
        sat_pos_earth = transfer_rv[0].to(u.AU)
        sat_pos_sun   = sat_pos_earth + earth_pos_sun

        # Get direction of Moon relative to the Spacecraft
        moon_pos_sat      = moon_pos_sun - sat_pos_sun
        moon_pos_sat_unit = moon_pos_sat/np.linalg.norm(moon_pos_sat)

        # Get direction of Star relative to the Spacecraft
        # Because the Star is really far, the direction is almost the same from the Sun or from the Spacecraft
        star1_pos_sat_unit = star1_pos_sun_unit
        star2_pos_sat_unit = star2_pos_sun_unit
        star3_pos_sat_unit = star3_pos_sun_unit

        # Compute the visible angle between the reference Star and the Moon
        alpha1[i] = np.arccos(np.dot(moon_pos_sat_unit,star1_pos_sat_unit))
        alpha2[i] = np.arccos(np.dot(moon_pos_sat_unit,star2_pos_sat_unit))
        alpha3[i] = np.arccos(np.dot(moon_pos_sat_unit,star3_pos_sat_unit))                                

plt.figure("CelesNav")
plt.plot(alpha1.to(u.deg),label='alpha 1')
plt.plot(alpha2.to(u.deg),label='alpha 2')
plt.plot(alpha3.to(u.deg),label='alpha 3')
plt.xlabel("Steps")
plt.ylabel("Angle (deg)")
plt.grid()
plt.legend()

plt.show()