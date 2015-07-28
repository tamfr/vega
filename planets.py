# -*- coding: utf-8 -*-
"""
Created on Mon Jul 27 21:20:51 2015

@author: Mott
"""
import numpy as np

muSun = 2.9591230378107664 * 10**(-4) * 149597870**3 / 86400**2 # Sun Standard Gravitational Parameter [km^3/s^2]

class planet(object):
    """Planet defines parameters of a planet:

    Attributes:
        e: Eccentricity.
        Ra: Aphelion [km].
        Rp: Perihelion [km]
        a: Semi-Major Axis [km].
        T: Orbital Period [s].
        r: Mean Radius [km].
        mu: Standard Graviational Parameter [km^3/s^2].
        g: Gravitational Acceleration [km/s^2].
    """
    def __init__(self):
        self.e = 0.093439031                              # Eccentricity
        self.Ra = 249234450.9                             # Aphelion [km]
        self.Rp = self.Ra * (1 - self.e) / (1 + self.e)   # Perihelion [km]
        self.a = self.Ra / (1 + self.e)                   # Semi-Major Axis [km]
        self.T = 2 * np.pi * (self.a**3 / muSun)**(1 / 2) # Orbital Period [s]
        self.r = 3389.92                                  # Mean Radius [km]
        self.mu = 42828.3                                 # Standard Graviational Parameter [km^3/s^2]
        self.g = 3.71 * 10**(-3)                          # Gravitational Acceleration [km/s^2]