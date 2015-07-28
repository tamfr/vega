# -*- coding: utf-8 -*-
"""
Created on Tue Jul 28 16:57:06 2015

@author: Mott
"""
import numpy as np

# Heliocentric Orbital System Constants (Source: NASA JPL Horizons Web Interface JDN: 2455562.5-2455927.5)

muSun = 2.9591230378107664*10**(-4)*149597870**3/86400**2; # Sun Standard Gravitational Parameter [km^3/s^2]
LonPerMars = 336.04084*np.pi/180; # Mars Longitude of Perihelion [radians] (http://nssdc.gsfc.nasa.gov/planetary/factsheet/earthfact.html)
LonPerEarth = 102.94719*np.pi/180; # Earth Longitude of Perihelion [radians] (http://nssdc.gsfc.nasa.gov/planetary/factsheet/earthfact.html)
Theta = LonPerMars-LonPerEarth; # Angle Between Earth's Perihelion and Mars's Perihelion [radians]
JDN0 = 2455565.2736; # Julian Day Number of Last Epoch (i.e. CE 2011 January 03 18:34:00 UT  Tuesday) 

##############################################################################
# Constant Orbital Elements for Earth

Earth = {
    'e' : 0.016741967,      # Eccentricity.
    'Rp': 152105805.7,      # Perihelion [km].
    'r' : 6378.136,         # Mean Radius [km].
    'mu': 398600.440,       # Standard Graviational Parameter [km^3/s^2].
    'g' : 9.81*10**(-3)     # Gravitational Acceleration [km/s^2].
}

##############################################################################
# Constant Orbital Elements for Mars

Mars = {
    'e' : 0.016741967,      # Eccentricity.
    'Rp': 152105805.7,      # Perihelion [km].
    'r' : 6378.136,         # Mean Radius [km].
    'mu': 398600.440,       # Standard Graviational Parameter [km^3/s^2].
    'g' : 9.81*10**(-3)     # Gravitational Acceleration [km/s^2].
}