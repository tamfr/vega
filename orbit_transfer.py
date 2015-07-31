# -*- coding: utf-8 -*-
"""
Created on Tue Jul 28 17:46:15 2015

@author: Mott
"""
# Important Remarks: 
# 1. J2000 Used to Establish All Constants
# 2. For Analysis, Earth's Perihelion is the Epoch (JDN 2455565.547 [i.e. CE 2011 January 04 01:07:40.8 UT  Tuesday])
# 3. Inclination and Orbital Precession are Neglected

from config import Earth, Mars
from planets import planet

Earth = planet(Earth)
Mars = planet(Mars)

Theta = Mars.lonPer-Earth.lonPer