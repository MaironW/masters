# Given:
# - Satellite position around the Earth (from orbit_generation)
# - Earth position around the Sun
# - Moon position around the Earth
# - Reference Star position (arbitrary vector)
# Get:
# - The visible angle between the reference Star and the Moon
# Optional:
# - Add Mars and Venus to the simulation
# - Assume other celestial bodies perturbation on the satellite dynamics

import numpy as np
import matplotlib.pyplot as plt

from astropy.time import Time
from astropy.coordinates import solar_system_ephemeris, get_body_barycentric
from poliastro.bodies import Sun, Earth, Moon
from poliastro.twobody import Orbit
from astropy.coordinates import CartesianRepresentation
import astropy.units as u

# Plot colors
c1 = (0.0,0.4,1.0)
c2 = (0.8,0.0,0.0)
c3 = (0.2,0.8,0.0)
c4 = 'm'

# Define direction of arbitrary Star relative to the Sun
# It must be inside the FOV of the Spacecraft to make sense
star_pos_sun_unit = [0.987, 0.101, 0.129] << u.one

# Time in Barycentric Dynamic Time
observation_time = Time("2025-04-27", scale="tdb")

# Set Spacecraft Orbit parameters
a    = 6731     << u.km
ecc  = 0.0001   << u.one
inc  = 51       << u.deg
raan = 13       << u.deg
argp = 42       << u.deg
nu   = 23.33    << u.deg # true anomaly

# Generate Keplerian Orbit
keplerian_orbit = Orbit.from_classical(Earth, a, ecc, inc, raan, argp, nu, epoch=observation_time)
kep_rv          = keplerian_orbit.rv()
sat_pos_earth   = kep_rv[0].to(u.AU)

with solar_system_ephemeris.set("builtin"):
    # Get Earth's position relative to the barycentric referential
    earth_pos_icrs = get_body_barycentric("Earth", observation_time)
    moon_pos_icrs  = get_body_barycentric("Moon",  observation_time)
    sun_pos_icrs   = get_body_barycentric("Sun",   observation_time)

# Convert positions to the Sun reference frame
earth_pos_sun = earth_pos_icrs - sun_pos_icrs
moon_pos_sun  = moon_pos_icrs  - sun_pos_icrs

# Convert to cartesian representation
earth_pos_sun = earth_pos_sun.represent_as(CartesianRepresentation).xyz.to(u.AU)
moon_pos_sun  = moon_pos_sun.represent_as(CartesianRepresentation).xyz.to(u.AU)

# Get position of Spacecraft relative to the Sun
sat_pos_sun = sat_pos_earth + earth_pos_sun

# Get direction of Moon relative to the Spacecraft
moon_pos_sat      = moon_pos_sun - sat_pos_sun
moon_pos_sat_unit = moon_pos_sat/np.linalg.norm(moon_pos_sat)

# Get direction of Star relative to the Spacecraft
# Because the Star is really far, the direction is almost the same from the Sun or from the Spacecraft
star_pos_sat_unit = star_pos_sun_unit

# Compute the visible angle between the reference Star and the Moon
alpha = np.arccos(np.dot(moon_pos_sat_unit,star_pos_sat_unit))

# Angle
print(alpha.to(u.deg))

# Plot
ax = plt.figure().add_subplot(projection='3d')
ax.plot(0,0,0,'s',color=c4,label='Spacecraft')
ax.plot([0,moon_pos_sat_unit[0]],[0,moon_pos_sat_unit[1]],[0,moon_pos_sat_unit[2]],'-',color=c1)
ax.plot(moon_pos_sat_unit[0],moon_pos_sat_unit[1],moon_pos_sat_unit[2],'o',color=c1,label='Moon')
ax.plot([0,star_pos_sat_unit[0]],[0,star_pos_sat_unit[1]],[0,star_pos_sat_unit[2]],'-',color=c2)
ax.plot(star_pos_sat_unit[0],star_pos_sat_unit[1],star_pos_sat_unit[2],'*',color=c2,label='Star')
ax.grid()
ax.legend()
plt.show()