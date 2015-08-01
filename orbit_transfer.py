# -*- coding: utf-8 -*-
"""
Created on Tue Jul 28 17:46:15 2015

@author: Mott
"""
# Important Remarks: 
# 1. J2000 Used to Establish All Constants
# 2. For Analysis, Earth's Perihelion is the Epoch (set in config settings)
# 3. Inclination and Orbital Precession are Neglected

from __future__ import division
from config import Earth, Mars, fE0, fM0, t_max, muSun, JDN0
from planets import planet
import orbital_tools as OT
import numpy as np

Earth = planet(Earth)
Mars = planet(Mars)

Theta = Mars.lonPer-Earth.lonPer # Angle between Mars perihelion and Earth perihelion.

T_E0 = OT.t_of_M_T( Earth.T, OT.M_of_E( Earth.e, OT.E_of_f( Earth.e, fE0 ) ) ) # Earth Period Advance at Epoch
T_M0 = OT.t_of_M_T( Mars.T,  OT.M_of_E( Mars.e,  OT.E_of_f( Mars.e,  fM0 ) ) ) # Mars Period Advance at Epoch

t = 0

while t < t_max:
    t = t + 86400 # Time increase [s] (86400 [s] is one day)
    
    # Planet positions over time.
    f_E = OT.f_of_E( Earth.e, OT.E_of_M( Earth.e, OT.M_of_t( Earth.T, T_E0 + t ) ) ) # Earth True anomaly given period advance from epoch
    f_M = OT.f_of_E( Mars.e,  OT.E_of_M( Mars.e,  OT.M_of_t( Mars.T,  T_M0 + t ) ) ) # Mars True anomaly given period advance from epoch
    
    # Hohmann Transfer Options
    
    f_MH = f_E - Theta + np.pi # Mars true anomaly on arrival for Hohmann transfer.
    R_MH = Mars.a*(1 - Mars.e**2)/(1 + Mars.e*np.cos(f_MH)) # Radial distance to Mars for Hohmann transfer.
    R_EH = Earth.a*(1 - Earth.e**2)/(1 + Earth.e*np.cos(f_E)) # Radial distance to Earth for Hohmann transfer.
    a_H = (R_MH + R_EH) / 2
    e_H = (R_MH - R_EH) / (R_MH + R_EH)
    T_H = np.pi*(a_H**3/muSun)**(1/2)
    
    f_MHA = OT.f_of_E( Mars.e,  OT.E_of_M( Mars.e,  OT.M_of_t( Mars.T,  T_M0 + t + T_H) ) ) # Mars true anomaly after desired Hohmann transfer time.

    
    if abs(f_MHA-f_MH) < 0.5*np.pi/180:
        print (JDN0*86400+t)/86400
        print T_H
        
    