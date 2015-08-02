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

def rising_transfer(f_POD, f_target, Theta, e_trans, a_trans, e_target, a_target, T_target, T_target0, e_POD, a_POD, T_POD, T_POD0, mu, return_path_option):
    lower = (a_target < a_POD) # Kind is either 0 for raising orbit and 1 for lowering orbit.
    
    f_TOA = OT.f_elliptic_transfer(f_POD, Theta, e_trans, a_trans, e_target, a_target, return_path_option, lower) # True anomaly of transfer on arrival.
        
    T_trans = OT.t_of_M_a(a_trans, mu, OT.M_of_E(e_trans, OT.E_of_f(e_trans, f_TOA))-np.cos(return_path_option*np.pi)*np.pi*lower) # Transfer time. 
        
    f_PTOA = OT.f_of_E( e_target,  OT.E_of_M( e_target,  OT.M_of_t( T_target,  T_target0 + t + T_trans) ) ) # Target planet true anomaly on arrival.
        
    if abs(f_PTOA-(f_TOA+f_POD-np.cos(lower*np.pi)*Theta-np.pi*lower)) < 0.5*np.pi/180:
        print 'Date of Fast Transfer [Julian Day Number]: \n' + str(JDN)  
                    
        f_PODOA = OT.f_of_E( e_POD,  OT.E_of_M( e_POD,  OT.M_of_t( T_POD,  T_POD0 + t + T_trans) ) ) # Earth true anomaly on Mars arrival for Hohmann transfer.
        
        print 'OT.transfer_plot('+str(f_POD) +', ' + str(f_PODOA) +', '+str(f_target)+', '+str(f_TOA+f_POD-np.cos(lower*np.pi)*Theta-np.pi*lower)+', '+str(Theta)+', '+ str(a_trans)+', '+str(e_trans)+',' +str(a_POD)+','+ str(e_POD)+', '+str(a_target)+', '+str(e_target)+', '+str(JDN)+')'

Earth = planet(Earth)
Mars = planet(Mars)

Theta = Mars.lonPer-Earth.lonPer # Angle between Mars perihelion and Earth perihelion.

T_E0 = OT.t_of_M_T( Earth.T, OT.M_of_E( Earth.e, OT.E_of_f( Earth.e, fE0 ) ) ) # Earth Period Advance at Epoch
T_M0 = OT.t_of_M_T( Mars.T,  OT.M_of_E( Mars.e,  OT.E_of_f( Mars.e,  fM0 ) ) ) # Mars Period Advance at Epoch

# Loop to solve for trajectories.

for t in xrange(0, t_max + step, step):
    JDN = (JDN0*86400+t)/86400 # Julian Day Number 
    #print JDN
    # Planet positions over time.
    f_E = OT.f_of_E( Earth.e, OT.E_of_M( Earth.e, OT.M_of_t( Earth.T, T_E0 + t ) ) ) # Earth True anomaly given period advance from epoch
    f_M = OT.f_of_E( Mars.e,  OT.E_of_M( Mars.e,  OT.M_of_t( Mars.T,  T_M0 + t ) ) ) # Mars True anomaly given period advance from epoch
    
    ##########################################
    ######## Hohmann Transfer Options ########
    
    f_MH = f_E - Theta + np.pi # Mars true anomaly on arrival for Hohmann transfer.
    R_MH = Mars.a*(1 - Mars.e**2)/(1 + Mars.e*np.cos(f_MH)) # Radial distance to Mars for Hohmann transfer.
    R_EH = Earth.a*(1 - Earth.e**2)/(1 + Earth.e*np.cos(f_E)) # Radial distance to Earth for Hohmann transfer.
    a_H = (R_MH + R_EH) / 2
    e_H = (R_MH - R_EH) / (R_MH + R_EH)
    T_H = np.pi*(a_H**3/muSun)**(1/2)
    
    f_MHA = OT.f_of_E( Mars.e,  OT.E_of_M( Mars.e,  OT.M_of_t( Mars.T,  T_M0 + t + T_H) ) ) # Mars true anomaly after desired Hohmann transfer time.

    
    if abs(f_MHA-f_MH) < 0.5*np.pi/180:
        
        f_EOMA = OT.f_of_E( Earth.e,  OT.E_of_M( Earth.e,  OT.M_of_t( Earth.T,  T_E0 + t + T_H) ) ) # Earth true anomaly on Mars arrival for Hohmann transfer.
        
        OT.transfer_plot(f_E, f_EOMA, f_M, f_MHA, Theta, a_H, e_H, Earth.a, Earth.e, Mars.a, Mars.e, JDN)
        
        print 'Date of Hohmann Departure [Julian Day Number]: \n' + str(JDN)        
        print 'Transfer Time [Days]:' + str(T_H/86400) + '\n'

    ############################################################   
    ######## All Outbound Trajectories to Target Planet ########
    
    e_trans = 0.999 # Highest eccentricity trnasfer to test for.
    
    e_step = (e_trans-e_H)/500 # Step size of Transfer Eccentricity
    
    while e_trans >= e_H:
        
        a_trans = R_EH/(1-e_trans)
        
        rising_transfer(f_E, f_M, Theta, e_trans, a_trans, Mars.e, Mars.a, Mars.T, T_M0, Earth.e, Earth.a, Earth.T, T_E0, muSun, 0)
        # Return Path options
        
        e_trans = e_trans - e_step
        
    ############################################################   
    ######## All Inbound Trajectories to Planet of Origin ########
        
#    f_EH = f_M + Theta - np.pi # Mars true anomaly on arrival for Hohmann transfer.
#    R_MH = Mars.a*(1 - Mars.e**2)/(1 + Mars.e*np.cos(f_M)) # Radial distance to Mars for Hohmann transfer.
#    R_EH = Earth.a*(1 - Earth.e**2)/(1 + Earth.e*np.cos(f_EH)) # Radial distance to Earth for Hohmann transfer.
#    a_H = (R_MH + R_EH) / 2
#    e_H = (R_MH - R_EH) / (R_MH + R_EH)
#    T_H = np.pi*(a_H**3/muSun)**(1/2)
#    
#    e_trans = 0.48
#    e_step = (e_trans-e_H)/500 # Step size of Transfer Eccentricity Earth Return
#    #nER = 0 # Initiate Options Per Day Counter
#    
#    while e_trans >= e_H:
#        a_trans = R_MH/(1+e_trans)
#        
#        #rising_transfer(f_M, f_E, Theta, e_trans, a_trans, Earth.e, Earth.a, Earth.T, T_E0, Mars.e, Mars.a, Mars.T, T_M0, muSun, 1)
#        Y = a_trans*(1-e_trans**2)/(Earth.a*(1-Earth.e**2))
#        f_T0 = np.pi/3#3*np.pi/2  # Initilize True Anomaly for Loop
#        f_T1 = -np.pi/2#-np.pi # Initilize True Anomaly for Loop
#        
#        while abs(f_T1-f_T0) > 0.0005: # Condition On Which to Run Loop
#            if f_T1 > np.pi:#2*np.pi:
#                f_T1 = f_T1-np.pi
#            elif f_T1 < 0:#np.pi:
#                f_T1 = f_T1+np.pi
#            else:
#                f_T0 = f_T1
#                f_MOA = f_T0+f_M+Theta-np.pi # np.cos(kind*np.pi) serves as a switch to make Theta negative and kind serves as an on/off switch for subtracting pi.
#                f_T1 = f_T0-(1+e_trans*np.cos(f_T0)-Y-Y*Earth.e*np.cos(f_MOA))/(Y*Earth.e*np.sin(f_MOA)-e_trans*np.sin(f_T0)) # Newton's Root Finding Method
#                
#        f_TOA = f_T1
#        
#        
#        T_trans = OT.t_of_M_a(a_trans, muSun, OT.M_of_E(e_trans, OT.E_of_f(e_trans, f_TOA))+np.pi) # Transfer time. 
#        
#        f_PTOA = OT.f_of_E( Earth.e,  OT.E_of_M( Earth.e,  OT.M_of_t( Earth.T,  T_E0 + t + T_trans) ) ) # Target planet true anomaly on arrival.
#        
#        if abs(f_PTOA-(f_TOA+f_M+Theta-np.pi)) < 0.5*np.pi/180:
#            print 'Date of return Transfer [Julian Day Number]: \n' + str(JDN)  
#            print 'Transfer time [days]: ' + str(T_trans/86400)        
#            f_PODOA = OT.f_of_E( Mars.e,  OT.E_of_M( Mars.e,  OT.M_of_t( Mars.T,  T_M0 + t + T_trans) ) ) # Earth true anomaly on Mars arrival for Hohmann transfer.
#        
#            print 'OT.plot_return('+str(f_M) +', ' + str(f_PODOA) +', '+str(f_E)+', '+str(f_TOA+f_M+Theta-np.pi)+', '+str(Theta)+', '+ str(a_trans)+', '+str(e_trans)+',' +str(Mars.a)+','+ str(Mars.e)+', '+str(Earth.a)+', '+str(Earth.e)+', '+str(JDN)+')'
#       
#       
#        e_trans = e_trans - e_step