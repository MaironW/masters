# Python script that converts orbital elements from state vector
# written by: Mairon S. Wolniewicz
# last update: 10/11/2024

import numpy as np
from numpy.linalg import norm
from numpy import dot
from numpy import cross
from numpy import arccos
from numpy import cos, sin
from numpy import dot
from numpy import pi

def elementsfromstate(pos,vel,mu):
    l,c = pos.shape

    sma  = [0]*l
    ecc  = [0]*l
    inc  = [0]*l
    argp = [0]*l
    raan = [0]*l
    f    = [0]*l

    for i in range(l):
        R = pos[i]
        V = vel[i]

        r      = norm(R)          # [m]      orbit radii
        v      = norm(V)          # [m/s]    orbit speed
        vr     = (dot(R,V))/r     # [m/s]    orbit radial speed
        H      = cross(R,V)       # [m^2/s]  angular momentum
        h      = norm(H)          # [m^2/s]  angular momentum norm
        inc[i] = arccos(H[2]/h)   # [rad]    orbit inclination
        N      = cross([0,0,1],H) # [   ]    node line direction
        n      = norm(N)          # [   ]    node line norm

        # if Nx/N>0, raan [I,IV]
        # if Nx/N<0, raan [II,III]
        raan[i] = arccos(N[0]/n) # [rad] right ascencion
        # if(N[1]<0):
            # raan[i] = 2*pi - raan[i]

        E = 1/mu * ((v**2 - mu/r)*R - r*vr*V) # eccentricity vector
        ecc[i] = norm(E) # eccentricity value

        # if N.E>0, argp [I,IV]
        # if N.E<0, argp [II,III]
        argp[i] = arccos(dot(N/n,E/ecc[i])) # [rad] argument of perigee
        if(E[2]<0):
            argp[i] = 2*pi - argp[i]

        # if E.R>0, f [I,IV]
        # if E.R<0, f [II,III]
        f[i] = arccos(dot(E/ecc[i],R/r)) # [rad] true anomaly
        if(vr<0):
            f[i] = 2*pi - f[i]

        q = h*h/mu * 1/(1+ecc[i]) # [km] perigee
        Q = h*h/mu * 1/(1-ecc[i]) # [km] apogee
        sma[i] = (q+Q)*0.5        # [km] semimajor axis
        T = 2*pi*mu**(-0.5) * sma[i]**(1.5) # [s] period

    return sma,ecc,inc,argp,raan,f

def statefromelements(a,e,I,argp,raan,f,mu=398600.0):
    h = (a*mu*(1-e**2))**0.5

    r = h**2/mu * 1/(1+e*cos(f)) * np.array([cos(f),sin(f),0])
    v = mu/h * np.array([-sin(f),e+cos(f),0])

    R1 = np.array([[cos(raan),sin(raan),0],
                   [-sin(raan),cos(raan),0],
                   [0,0,1]])

    R2 = np.array([[1,0,0],
                   [0,cos(I),sin(I)],
                   [0,-sin(I),cos(I)]])

    R3 = np.array([[cos(argp),sin(argp),0],
                   [-sin(argp),cos(argp),0],
                   [0,0,1]])

    r = dot(dot(dot(r,R3),R2),R1)
    v = dot(dot(dot(v,R3),R2),R1)

    return r,v