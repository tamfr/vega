# -*- coding: utf-8 -*-
"""
Created on Mon Jul 27 21:20:51 2015

@author: Mott
"""
from __future__ import division
import numpy as np
import orbital_tools as OT

# Sun Standard Gravitational Parameter [km^3/s^2]

muSun = 2.9591230378107664 * 10**(-4) * 149597870**3 / 86400**2

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
    def __init__(self, planet):
        self.e = planet['e']                                # Eccentricity
        self.Ra = planet['Ra']                              # Aphelion [km]
        self.Rp = self.Ra * (1 - self.e) / (1 + self.e)     # Perihelion [km]
        self.a = self.Ra / (1 + self.e)                     # Semi-Major Axis [km]
        self.T = 2 * np.pi * (self.a**3 / muSun)**(1 / 2)   # Orbital Period [s]
        self.r = planet['r']                                # Mean Radius [km]
        self.mu = planet['mu']                              # Standard Graviational Parameter [km^3/s^2]
        self.g = planet['g']                                # Gravitational Acceleration [km/s^2]
        self.lonPer = planet['lonPer']                      # Longitude of Perihelion [radians] 
        self.orbit = OT.orbit(self.a, self.e)
    
    def M(self, t):    
        M = OT.M_of_t( self.T, t )    
        return M
        
    def E(self, t):    
        M = OT.M_of_t( self.T, t )
        E = OT.E_of_M( self.e, M )
        return E
        
    def f(self, t):
        M = OT.M_of_t( self.T, t )
        E = OT.E_of_M( self.e, M )
        f = OT.f_of_E( self.e, E )
        return f
