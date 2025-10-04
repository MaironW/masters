# Load SPICE Kernels

import spiceypy

_KERNELS_LOADED = False

def load_kernels():
    global _KERNELS_LOADED
    if _KERNELS_LOADED:
        return
    kernel_dir = "Utils/kernels/"
    spiceypy.furnsh(kernel_dir + "naif0012.tls") # Leap seconds
    spiceypy.furnsh(kernel_dir + "de440s.bsp")   # Planetary ephemerides
    spiceypy.furnsh(kernel_dir + "mar097.bsp")   # Mars body + moons
    spiceypy.furnsh(kernel_dir + "pck00010.tcp") # Planet orientation

    _KERNELS_LOADED = True

def clear_kernels():
    spiceypy.kclear()

# Convert a time string to J2000 epoch in seconds
def get_time(time: str):
    load_kernels()
    return spiceypy.str2et(time)

# Give the position and velocity of the body defined in the SSB
def get_state(body: str, et: float):
    load_kernels()
    # body:   ["SUN","EARTH","MARS", etc.]
    # et:     Observer epoch in seconds past J2000
    # ref:    Reference frame of output state vector
    # abcorr: Aberration correction flag
    # obs:    Observing body frame
    # state:  pos [km], vel [km/s]
    # lt:     One way light time between observer and target [s]
    state, lt = spiceypy.spkezr(body, et, "J2000", "NONE", "SOLAR SYSTEM BARYCENTER")
    return state[:3], state[3:] # pos, vel

# Give the orientation of the BCBF (Body-Fixed, Body-Centered) referencial with respect to the SSB
def get_orientation(body: str, et: float):
    load_kernels()
    frame_name = f"IAU_{body.upper()}"
    # Return the matrix that transforms position vectors from one frame to another at a specified epoch
    r = spiceypy.pxform("J2000", frame_name, et)
    # Convert matrix to quaternion
    return spiceypy.m2q(r)
