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
from config import Earth, Mars, fE0, fM0, t_max, muSun, JDN0, step
from planets import planet
import orbital_tools as OT
import numpy as np
import matplotlib.pyplot as plt

Earth = planet(Earth)
Mars = planet(Mars)

Theta = Mars.lonPer-Earth.lonPer # Angle between Mars perihelion and Earth perihelion.

T_E0 = OT.t_of_M_T( Earth.T, OT.M_of_E( Earth.e, OT.E_of_f( Earth.e, fE0 ) ) ) # Earth Period Advance at Epoch
T_M0 = OT.t_of_M_T( Mars.T,  OT.M_of_E( Mars.e,  OT.E_of_f( Mars.e,  fM0 ) ) ) # Mars Period Advance at Epoch

# Loop to solve for trajectories.

for t in xrange(0, t_max + step, step):
    
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
        
        f_EOMA = OT.f_of_E( Earth.e,  OT.E_of_M( Earth.e,  OT.M_of_t( Earth.T,  T_E0 + t + T_H) ) ) # Earth true anomaly on Mars arrival for Hohmann transfer.
        
        
        JDN = (JDN0*86400+t)/86400
        print T_H
        
        fig = plt.figure() # initialize figure
        ax = fig.add_subplot(111) # name of the plot
        f = np.arange(0,2*np.pi,2*np.pi/999)
        ax.set_aspect('equal')
        
        R_T = a_H*(1 - e_H**2)/(1 + e_H*np.cos(f))
        R_E = Earth.a*(1 - Earth.e**2)/(1 + Earth.e*np.cos(f))
        R_M = Mars.a*(1 - Mars.e**2)/(1 + Mars.e*np.cos(f))
            
        ax.plot(R_T*np.cos(f+f_E), R_T*np.sin(f+f_E), color='g', lw=1)     # Plots the transfer's orbit in green.
        ax.plot(R_E*np.cos(f), R_E*np.sin(f), color='b', lw=1)             # Plots Earth's orbit in blue.
        ax.plot(R_M*np.cos(f+Theta), R_M*np.sin(f+Theta), color='r', lw=1) # Plots Mars's orbit in red.
        ax.plot(OT.R(Earth.a, Earth.e, f_E)*np.cos(f_E), OT.R(Earth.a, Earth.e, f_E)*np.sin(f_E), color='b', marker='o')              # Plots Earth's position upon TMI.
        ax.plot(OT.R(Mars.a, Mars.e, f_M)*np.cos(f_M+Theta), OT.R(Mars.a, Mars.e, f_M)*np.sin(f_M+Theta), color='r', marker='o')      # Plots Mars's position upon TMI.       
        ax.plot(OT.R(Earth.a, Earth.e, f_EOMA)*np.cos(f_EOMA), OT.R(Earth.a, Earth.e, f_EOMA)*np.sin(f_EOMA), color='b', marker='o', fillstyle='none')  # Plots Earth's position upon MOI.        
        ax.plot(OT.R(Mars.a, Mars.e, f_MH)*np.cos(f_MH+Theta), OT.R(Mars.a, Mars.e, f_MH)*np.sin(f_MH+Theta), color='g', marker='o', fillstyle='none')  # Plots Mars's position upon MOI.
        ax.plot(0, 0, color='k', marker='x') # Plots a black "x" to indicate the Sun's location.

        plt.savefig("result" + str(JDN) + ".eps", format="eps")        
        plt.show()
        