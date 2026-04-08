# Miscelaneous utility functions used throughout the code

import numpy as np

# Convert from Latitude and Longitude to XYZ
def latlon2dir(lat, lon):
    lat = np.deg2rad(lat)
    lon = np.deg2rad(lon)

    clat, slat = np.cos(lat), np.sin(lat)
    clon, slon = np.cos(lon), np.sin(lon)

    dir_vec = np.array([clat*clon, clat*slon, slat])

    return dir_vec

# Convert vector from Galactic frame to SSB frame
def GALtoSSB(vec_GAL):
    # Rotation matrix from Galactic to Equatorial plane J2000
    R = np.array([
        [-0.0548755604, -0.8734370902, -0.4838350155],
        [+0.4941094279, -0.4448296300, +0.7469822445],
        [-0.8676661490, -0.1980763734, +0.4559837762]
    ])
    vec_SSB = (R @ vec_GAL).T
    return vec_SSB

# Generate a gnomonic projection of a vector
def gnomonic_projection(vec):
    x = vec[:,0]
    y = vec[:,1]
    z = vec[:,2]
    mask = z > 0
    u = np.full_like(z, np.nan)
    v = np.full_like(z, np.nan)
    u[mask] = x[mask] / z[mask]
    v[mask] = y[mask] / z[mask]
    return u, v

# Return the boundary of the Gnomonic projection
def gnomonic_boundary(field_of_view):
    phi = np.linspace(0, 2*np.pi, 500)
    r = np.tan(field_of_view)
    u = r * np.cos(phi)
    v = r * np.sin(phi)
    return u, v

# Generate Aitoff projection of a vector, in degrees
def aitoff_projection(vec):
    x = vec[:,0]
    y = vec[:,1]
    z = vec[:,2]

    ra   = np.arctan2(y, x)
    decl = np.arcsin(z)

    theta = np.arccos(np.cos(decl)*np.cos(ra/2))

    sinc_theta = np.where(
        np.abs(theta) < 1e-12,
        1.0,
        np.sin(theta)/theta
    )

    u = 2*np.cos(decl)*np.sin(ra/2)/sinc_theta
    v = np.sin(decl)/sinc_theta
    u = np.rad2deg(u)
    v = np.rad2deg(v)
    return u, v

# Return the boundary of the Aitoff projection, in degrees
def aitoff_boundary():
    t  = np.linspace(0, 2*np.pi, 500)
    x = np.pi*np.cos(t)
    y = (np.pi/2)*np.sin(t)
    x = np.rad2deg(x)
    y = np.rad2deg(y)
    return x, y

# Return true if array is entirely composed by nans
def is_nan(array):
    return np.isnan(array).all()
