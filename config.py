# -*- coding: utf-8 -*-
"""
Created on Tue Jul 28 16:57:06 2015

@author: Mott
"""
import numpy as np

# Heliocentric Orbital System Constants (Source: NASA JPL Horizons Web Interface JDN: 2455562.5-2455927.5)

muSun = 2.9591230378107664*10**(-4)*149597870**3/86400**2; # Sun Standard Gravitational Parameter [km^3/s^2]

##############################################################################
# Constant Orbital Elements for Earth

Earth = {
    'e'     : 0.016741967,          # Eccentricity.
    'Ra'    : 152105805.7,          # Aphelion [km].
    'r'     : 6378.136,             # Mean Radius [km].
    'mu'    : 398600.440,           # Standard Graviational Parameter [km^3/s^2].
    'g'     : 9.81*10**(-3),        # Gravitational Acceleration [km/s^2].
    'lonPer': 102.94719*np.pi/180   # Longitude of Perihelion [radians] (http://nssdc.gsfc.nasa.gov/planetary/factsheet/earthfact.html)
}

##############################################################################
# Constant Orbital Elements for Mars

Mars = {
    'e'     : 0.093439031,          # Eccentricity.
    'Ra'    : 249234450.9,          # Aphelion [km].
    'r'     : 3389.92,              # Mean Radius [km].
    'mu'    : 42828.3,              # Standard Graviational Parameter [km^3/s^2].
    'g'     : 3.71*10**(-3),        # Gravitational Acceleration [km/s^2].
    'lonPer': 336.04084*np.pi/180   # Longitude of Perihelion [radians] (http://nssdc.gsfc.nasa.gov/planetary/factsheet/earthfact.html)
}

# Theta = LonPerMars-LonPerEarth # Angle Between Earth's Perihelion and Mars's Perihelion [radians]

# Source for perihelion data: http://www.astropixels.com/ephemeris/perap2001.html and NASA Horizons.

# fE0 = (4.54323446*10**(-4))*np.pi/180; fM0 = 319.4234417*pi/180; JDN0 = 2455565.2736 # Date of Epoch: 03 Jan 2011 18:34:00 UT
# fE0 = (7.8331757*10**(-1))*np.pi/180; fM0 = 199.18562*pi/180; JDN0 = 2457391.2736 # Date of Epoch: 03 Jan 2016 18:34:00 UT
#fE0 = 3.5999934*10**2*np.pi/180; fM0 = 9.0419941*10**1*np.pi/180; JDN0 = 2459217.0770 # Date of Epoch: 02 Jan 2021 13:51:00 UT
fE0 = 3.5999922*10**2*np.pi/180; fM0 = 3.093210195*10**2*np.pi/180; JDN0 = 2461044.2194 # Date of Epoch: 03 Jan 2026 17:16:00 UT
# fE0 = (3.5988321*10**(2))*np.pi/180; fM0 = 192.14974*np.pi/180; JDN0 = 2462871.2736 # Date of Epoch: 04 Jan 2031 18:34:00 UT
 
t_max = 86400*365*3 # Length of search window [s]
step = 86400 # Step size [s] (86400 [s] is one day)

round_trip_max_days = 720 # Maximum round trip time [days]
min_stay_time = 15 # Minimum stay time on target planet [days]

# Mission requirements
alt_park_Earth = 330 # Parking orbit altitude [km]

alt_park_Mars = 150 # Parking orbit altitude [km]