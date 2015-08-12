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
from config import Earth, Mars, fE0, fM0, t_max, muSun, JDN0, step, round_trip_max_days, min_stay_time, font
from planets import planet
import orbital_tools as OT
import numpy as np

def elliptical_transfer(
    f_POD, 
    f_target, 
    Theta, 
    e_trans, 
    a_trans, 
    target,  
    T_target0, 
    POD, 
    T_POD0, 
    mu, 
    return_path_option = 0,
    ):

    """Elliptical transfer.
        Requires: 
            f_POD: true anomaly of planet of departure at time of transfer;
            f_target: true anomaly of target planet at time of transfer;           
            Theta: difference between target orbit and orbit of departure longitude of perihelion;            
            e_trans: transfer orbit eccentricity;
            a_trans: transfer orbit semi-major axis;
            target: target planet (instance);
            T_target0: target orbit period advance at epoch;
            POD: planet of departure (instance);
            T_POD0: planet of departure period advnace at epoch;
            mu: standard gravitational parameter of central body;
            return_path_option: set 1 for yes and 0 for no.
    """
    
    lower = (target.a < POD.a) # Kind is either 0 for raising orbit and 1 for lowering orbit.
    
    f_TOA = OT.f_elliptic_transfer(f_POD, Theta, e_trans, a_trans, target.e, target.a, return_path_option, lower) # True anomaly of transfer on arrival.
        
    T_trans = OT.t_of_M_a(a_trans, mu, OT.M_of_E(e_trans, OT.E_of_f(e_trans, f_TOA))-np.cos(return_path_option*np.pi)*np.pi*lower) # Transfer time. 
        
    f_PTOA = OT.f_of_E( target.e,  OT.E_of_M( target.e,  OT.M_of_t( target.T,  T_target0 + t + T_trans) ) ) # Target planet true anomaly on arrival.
        
    if abs(f_PTOA-(f_TOA+f_POD-np.cos(lower*np.pi)*Theta-np.pi*lower)) < 0.5*np.pi/180:
        f_PODOA = OT.f_of_E( POD.e,  OT.E_of_M( POD.e,  OT.M_of_t( POD.T,  T_POD0 + t + T_trans) ) ) # Earth true anomaly on Mars arrival for Hohmann transfer.
        if lower == 0: 
            #print 'Date of outbound transfer [Julian Day Number]: \n' + str(JDN)  
            output = 'OT.transfer_plot('+str(f_POD) +', ' + str(f_PODOA) +', '+str(f_target)+', '+str(f_TOA+f_POD-np.cos(lower*np.pi)*Theta-np.pi*lower)+', '+str(Theta)+', '+ str(a_trans)+', '+str(e_trans)+',' +str(POD.a)+','+ str(POD.e)+', '+str(target.a)+', '+str(target.e)+', '+str(JDN)+')'
        else:
            #print 'Date of return transfer [Julian Day Number]: \n' + str(JDN)            
            output = 'OT.plot_return('+str(f_M) +', ' + str(f_PODOA) +', '+str(f_E)+', '+str(f_TOA+f_M+Theta-np.pi)+', '+str(Theta)+', '+ str(a_trans)+', '+str(e_trans)+',' +str(Mars.a)+','+ str(Mars.e)+', '+str(Earth.a)+', '+str(Earth.e)+', '+str(JDN)+')'
        
        return np.array([JDN, JDN+T_trans/86400, T_trans, a_trans, e_trans, POD, target, f_E, f_M, f_TOA, f_PTOA, f_PODOA, return_path_option, lower])

Earth = planet(Earth)
Mars = planet(Mars)

Theta = Mars.lonPer-Earth.lonPer # Angle between Mars perihelion and Earth perihelion.

 # Earth Period Advance at Epoch

E = OT.E_of_f( Earth.e, fE0 )
M = OT.M_of_E( Earth.e, E )
T_E0 = OT.t_of_M_T( Earth.T, M )

# Mars Period Advance at Epoch

E = OT.E_of_f( Mars.e,  fM0 )
M = OT.M_of_E( Mars.e,  E) 
T_M0 = OT.t_of_M_T( Mars.T,  M ) 

out_transfers = np.zeros([14])
return_transfers = np.zeros([14])

# Loop to solve for trajectories.

for t in xrange(0, t_max + step, step):
    JDN = (JDN0*86400+t)/86400 # Julian Day Number 
    #print JDN
    # Planet positions over time.
    f_E = Earth.f(T_E0 + t) # Earth True anomaly given period advance from epoch
    f_M = Mars.f(T_M0 +t) # Mars True anomaly given period advance from epoch
    
    ##########################################
    ######## Hohmann Transfer Options ########
    
    f_MH = f_E - Theta + np.pi                                 # Mars true anomaly on arrival for Hohmann transfer.
    R_MH = Mars.a*(1 - Mars.e**2)/(1 + Mars.e*np.cos(f_MH))    # Radial distance to Mars for Hohmann transfer.
    R_EH = Earth.a*(1 - Earth.e**2)/(1 + Earth.e*np.cos(f_E))  # Radial distance to Earth for Hohmann transfer.
    a_H = (R_MH + R_EH) / 2
    e_H = (R_MH - R_EH) / (R_MH + R_EH)
    T_H = np.pi*(a_H**3 / muSun)**(1/2)
    
    f_MHA = OT.f_of_E( Mars.e,  OT.E_of_M( Mars.e,  OT.M_of_t( Mars.T,  T_M0 + t + T_H) ) ) # Mars true anomaly after desired Hohmann transfer time.

    
    if abs(f_MHA-f_MH) < 0.5*np.pi/180:
        
        M = OT.M_of_t( Earth.T,  T_E0 + t + T_H)        
        E = OT.E_of_M( Earth.e,  M )
        f_EOMA = OT.f_of_E( Earth.e,  E ) # Earth true anomaly on Mars arrival for Hohmann transfer.
        
        OT.transfer_plot(f_E, f_EOMA, f_M, f_MHA, Theta, a_H, e_H, Earth, Mars, JDN, font)
        
        print 'Date of Hohmann Departure [Julian Day Number]: \n' + str(JDN)        
        print 'Transfer Time [Days]:' + str(T_H/86400) + '\n'

    ############################################################   
    ######## All Outbound Trajectories to Target Planet ########
    
    e_trans = 0.999 # Highest eccentricity trnasfer to test for.
    
    e_step = (e_trans-e_H)/500 # Step size of Transfer Eccentricity
    
    while e_trans >= e_H:
        
        a_trans = R_EH/(1-e_trans)
        
        transfer = elliptical_transfer(f_E, f_M, Theta, e_trans, a_trans, Mars, T_M0, Earth, T_E0, muSun, 0)
        if str(transfer) != 'None':        
            out_transfers = np.vstack((out_transfers, transfer))        
        
        # Return path options (i.e. The transfer orbit picks up the target planet on second pass of the target's orbit. )
        transfer = elliptical_transfer(f_E, f_M, Theta, e_trans, a_trans, Mars, T_M0, Earth, T_E0, muSun, 1)
        if str(transfer) != 'None':        
            out_transfers = np.vstack((out_transfers, transfer))        
        
        e_trans = e_trans - e_step
        
    ############################################################   
    ######## All Inbound Trajectories to Planet of Origin ########
        
    f_EH = f_M + Theta - np.pi # Mars true anomaly on arrival for Hohmann transfer.
    R_MH = Mars.a*(1 - Mars.e**2)/(1 + Mars.e*np.cos(f_M)) # Radial distance to Mars for Hohmann transfer.
    R_EH = Earth.a*(1 - Earth.e**2)/(1 + Earth.e*np.cos(f_EH)) # Radial distance to Earth for Hohmann transfer.
    a_H = (R_MH + R_EH) / 2
    e_H = (R_MH - R_EH) / (R_MH + R_EH)
    T_H = np.pi*(a_H**3/muSun)**(1/2)
    
    e_trans = 0.48
    e_step = (e_trans-e_H)/500 # Step size of Transfer Eccentricity Earth Return
#    #nER = 0 # Initiate Options Per Day Counter
    
    while e_trans >= e_H:
        a_trans = R_MH/(1+e_trans)
        
        transfer = elliptical_transfer(f_M, f_E, Theta, e_trans, a_trans, Earth, T_E0, Mars, T_M0, muSun, 1)
        if str(transfer) != 'None':        
            return_transfers = np.vstack((return_transfers, transfer))
        
        e_trans = e_trans - e_step

mission_profiles=np.zeros([28])
for transfer in out_transfers:
    for return_transfer in return_transfers:
        if return_transfer[1] - transfer[0] < round_trip_max_days and  return_transfer[0] - transfer[1] > min_stay_time:
            mission_profiles = np.vstack((mission_profiles, np.hstack((transfer, return_transfer))))
        