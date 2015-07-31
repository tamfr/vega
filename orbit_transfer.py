# -*- coding: utf-8 -*-
"""
Created on Tue Jul 28 17:46:15 2015

@author: Mott
"""
# Important Remarks: 
# 1. J2000 Used to Establish All Constants
# 2. For Analysis, Earth's Perihelion is the Epoch (JDN 2455565.547 [i.e. CE 2011 January 04 01:07:40.8 UT  Tuesday])
# 3. Inclination and Orbital Precession are Neglected

from config import Earth, Mars, fE0, muSun
from planets import planet
import orbital_tools as OT
import numpy as np

Earth = planet(Earth)
Mars = planet(Mars)

Theta = Mars.lonPer-Earth.lonPer




ME0 = OT.M_of_E(Earth.e,OT.E_of_f(Earth.e,fE0))
TE0 = ME0/(2*np.pi/Earth.T) # Earth Period Advance at Epoch
f = OT.f_of_E(Earth.e,OT.E_of_M(Earth.e,OT.M_of_t(Earth.T,TE0 + 365*86400)))
