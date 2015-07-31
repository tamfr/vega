# -*- coding: utf-8 -*-
"""
Created on Tue Jul 28 17:46:15 2015

@author: Mott
"""
# Important Remarks: 
# 1. J2000 Used to Establish All Constants
# 2. For Analysis, Earth's Perihelion is the Epoch (set in config settings)
# 3. Inclination and Orbital Precession are Neglected

from config import Earth, Mars, fE0, fM0, t_max, muSun
from planets import planet
import orbital_tools as OT
import numpy as np

Earth = planet(Earth)
Mars = planet(Mars)

Theta = Mars.lonPer-Earth.lonPer

TE0 = OT.t_of_M_T(Earth.T,OT.M_of_E(Earth.e,OT.E_of_f(Earth.e,fE0))) # Earth Period Advance at Epoch
TM0 = OT.t_of_M_T(Mars.T,OT.M_of_E(Mars.e,OT.E_of_f(Mars.e,fM0))) # Mars Period Advance at Epoch

t = 0

while t < t_max:
    t = t + 84600
    
    fE = OT.f_of_E(Earth.e,OT.E_of_M(Earth.e,OT.M_of_t(Earth.T,TE0 + t))) # True anomaly given period advance from epoch
    fM = OT.f_of_E(Mars.e,OT.E_of_M(Mars.e,OT.M_of_t(Mars.T,TM0 + t)))
    #print(fE)