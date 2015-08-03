# -*- coding: utf-8 -*-
"""
Created on Tue Jul 28 22:09:52 2015

@author: Mott
"""
from __future__ import division
import numpy as np
import matplotlib.pyplot as plt

def R(a,e,f):
    return a*(1 - e**2)/(1 + e*np.cos(f))

####################### Orbital elements given time #######################

def M_of_t(T,t): 
    """Mean anomaly as fucntion of time.
        Requires period.
    """    
    return (2*np.pi/T)*((t)-int((t)/T)*T)    

def E_of_M(e,M): 
    """Eccentric anomaly as fucntion of mean anomaly.
        Requires eccentritcity.
    """
    count = 0
    # Iteration Loop to Find Eccentric Anomaly
    
    E0 = -np.pi # Initilize Eccentric Anomaly for Loop
    E1 = M+(e+np.cos(M))*np.sin(M) # Initilize Eccentric Anomaly for Loop
    
    while abs(E1-E0) > 0.0005: # Condition On Which to Run Loop
        count += 1        
        if E1 > 2*np.pi:
            E1 = E1 - 2*np.pi
        elif E1 < 0:
            E1 = E1 + np.pi
        else:
            E0 = E1
            E1 = E0 - (E0 - e*np.sin(E0) - M) / (1 - e*np.cos(E0)) # Newton's Root Finding Method
    return E1

def E_of_M_old(e,M):    
    """Eccentric anomaly as fucntion of mean anomaly.
        Requires eccentritcity.
    """
    count = 0
    # Iteration Loop to Find Eccentric Anomaly
    
    E0 = np.pi/3 # Initilize Eccentric Anomaly for Loop
    E1 = -np.pi/2 # Initilize Eccentric Anomaly for Loop
    
    while abs(E1-E0) > 0.0005: # Condition On Which to Run Loop
        count += 1        
        if E1 > 2*np.pi:
            E1 = E1 - 2*np.pi
        elif E1 < 0:
            E1 = E1 + np.pi
        else:
            E0 = E1
            E1 = E0 - (E0 - e*np.sin(E0) - M) / (1 - e*np.cos(E0)) # Newton's Root Finding Method
    return (E1,count)
    
def f_of_E(e,E):
    """True anomaly as fucntion of eccentric anomaly.
        Requires eccentritcity.
    """
    f = 2*np.arctan(((1 + e) / (1 - e))**(1/2)*np.tan(E/2)) # True Anomaly
    
    # Test to Make Sure True Anomaly is Between 0 and 360 Degrees
    if f < 0: 
        f = f + 2*np.pi
    elif f > 2*np.pi:
        f = f - 2*np.pi
    
    return f
    
####################### Orbital elements given true anomaly #######################

def E_of_f(e,f):
    """Eccentric anomaly as fucntion of true anomaly.
        Requires eccentritcity.
    """    
    E = 2*np.arctan(((1-e)/(1+e))**(1/2)*np.tan(f/2)) # Determines Transfer Orbit Eccentric Anomaly    
    # Test to Make Sure Eccentric Anomaly is Greater than 0.
    if E < 0:
        E = E + 2*np.pi
        
    return E

def M_of_E(e,E):
    """Mean anomaly as fucntion of eccentric anomaly.
        Requires eccentritcity.
    """ 
    return E-e*np.sin(E)
    
def t_of_M_a(a,mu,M):
    """Mean anomaly as fucntion of eccentric anomaly.
        Requires semi-major axis and standard gravitational parameter of the central body.
    """ 
    return M*(a**3/mu)**(1/2)

def t_of_M_T(T,M):
    """Mean anomaly as fucntion of eccentric anomaly.
        Requires semi-major axis and standard gravitational parameter of the central body.
    """ 
    return M*T/(2*np.pi)


####################### Non-Hohmann Elliptical Transfers #######################

def f_elliptic_transfer(f, Theta, e_trans, a_trans, e_target, a_target, return_path_option, lower):
    """Elliptical transfer true anomaly.
        Requires: 
            f: True anomaly of planet of departure at time of transfer;
            Theta: Difference between target orbit and orbit of departure longitude of perihelion;            
            e_trans: transfer orbit eccentricity;
            a_trans: transfer orbit semi-major axis;
            e_target: target orbit eccentricity;
            a_target: target orbit semi-major axis;
            mu: standard gravitational parameter of the central body.
    """       
    low = (not return_path_option and lower)
    high = (return_path_option and not lower)    
    
    
    Y = a_trans*(1-e_trans**2)/(a_target*(1-e_target**2))
    if lower == 0:
        f_T0 = np.pi/3*(not return_path_option) + 3*np.pi/2*(return_path_option)  # Initilize True Anomaly for Loop
        f_T1 = -np.pi/2*(not return_path_option)  + 2*np.pi/3*(return_path_option) # Initilize True Anomaly for Loop
    
    if lower == 1:
        f_T0 = 3*np.pi/2*(not return_path_option) + np.pi/3*( return_path_option)  # Initilize True Anomaly for Loop
        f_T1 = -np.pi*(not return_path_option)  + -np.pi/3*(return_path_option)# Initilize True Anomaly for Loop

    while abs(f_T1-f_T0) > 0.0005: # Condition On Which to Run Loop
        if f_T1 > np.pi*(((not return_path_option and not lower) + (return_path_option and lower)) + 2*(high+low)):
            f_T1 = f_T1-np.pi
        elif f_T1 < np.pi*(high + low):
            f_T1 = f_T1+np.pi
        else:
            f_T0 = f_T1
            f_MOA = f_T0+f-np.cos(lower*np.pi)*Theta-np.pi*lower # np.cos(lower*np.pi) serves as a switch to make Theta negative and kind serves as an on/off switch for subtracting pi.
            f_T1 = f_T0-(1+e_trans*np.cos(f_T0)-Y-Y*e_target*np.cos(f_MOA))/(Y*e_target*np.sin(f_MOA)-e_trans*np.sin(f_T0)) # Newton's Root Finding Method
    return f_T1
    
#def f_slow_transfer(f, Theta, e_trans, a_trans, e_target, a_target):
#    """Elliptical fast transfer.
#        Requires: 
#            f: True anomaly of planet of departure at time of transfer;
#            Theta: Difference between target orbit and orbit of departure longitude of perihelion;            
#            e_trans: transfer orbit eccentricity;
#            a_trans: transfer orbit semi-major axis;
#            e_target: target orbit eccentricity;
#            a_target: target orbit semi-major axis.
#    """       
#    Y = a_trans*(1-e_trans**2)/(a_target*(1-e_target**2))
#    f_T0 = 3*np.pi/2  # Initilize True Anomaly for Loop
#    f_T1 = 2*np.pi/3 # Initilize True Anomaly for Loop
#    
#    while abs(f_T1-f_T0) > 0.0005: # Condition On Which to Run Loop
#     
#        if f_T1 > 2*np.pi:
#            f_T1 = f_T1-np.pi
#        elif f_T1 < np.pi:
#            f_T1 = f_T1+np.pi
#        else:
#            f_T0 = f_T1
#            f_MOA = f_T0+f-Theta
#            f_T1 = f_T0-(1+e_trans*np.cos(f_T0)-Y-Y*e_target*np.cos(f_MOA))/(Y*e_target*np.sin(f_MOA)-e_trans*np.sin(f_T0)) # Newton's Root Finding Method
#        
#    return f_T1

####################### Delta V calculations ####################### 
   
def delta_V_approx_TPI(e_trans, a_trans, a_POD, R_POD, mu_central_body):
    """Approximate Delta V for trans-planet injection.
        Requires: 
            e_trans: transfer orbit eccentricity;
            a_trans: transfer orbit semi-major axis;
            a_POD: Planet of departure semi-major axis;
            R_POD: Planet of departure radial distance at time of departure.
    """
    return  (((1+e_trans)/(1-e_trans))*mu_central_body/a_trans)**(1/2)-(mu_central_body*(2/R_POD-1/a_POD))**(1/2) # Delta V for trans-planet injection [km/s].  
          
def delta_V_actual_TPI(alt_park, mu_POD, r_POD, delta_v_approx_TPI):
    """Actual Delta V for trans-planet injection.
         Requires: 
            alt_park: parking orbit altitude of planet of departure;
            mu_POD: standard gravitational parameter of planet of departure;
            r_POD: radius of planet of departure;
            delta_v_approx_TPI: Approximate delta V for trans planet injection.
    """     
    PerHyper = alt_park + r_POD # Periapsis Distance From Center of Earth [km]
    a_h = mu_POD/delta_v_approx_TPI**2 # Determines Hyperbolic Semi-Major Axis [km]
    e_h = PerHyper/a_h + 1 # Determines Hyperbolic Eccentricity
    return ((e_h + 1)/(e_h - 1))^(1/2)*delta_v_approx_TPI-(mu_POD/PerHyper)^(1/2) # Actual Delta V for trans-planet injection [km/s]
 
####################### Plotting orbit transfers #######################
   
def transfer_plot(f_POD, f_PODOA, f_target, f_targetOA, Theta, a_trans, e_trans, a_POD, e_POD, a_target, e_target, JDN):
        fig = plt.figure() # initialize figure
        ax = fig.add_subplot(111) # name of the plot
        f = np.arange(0,2*np.pi,2*np.pi/999)
        ax.set_aspect('equal')
        
        R_T = a_trans*(1 - e_trans**2)/(1 + e_trans*np.cos(f))
        R_E = a_POD*(1 - e_POD**2)/(1 + e_POD*np.cos(f))
        R_M = a_target*(1 - e_target**2)/(1 + e_target*np.cos(f))
            
        ax.plot(R_T*np.cos(f+f_POD), R_T*np.sin(f+f_POD), color='g', lw=1)     # Plots the transfer's orbit in green.
        ax.plot(R_E*np.cos(f), R_E*np.sin(f), color='b', lw=1)             # Plots Earth's orbit in blue.
        ax.plot(R_M*np.cos(f+Theta), R_M*np.sin(f+Theta), color='r', lw=1) # Plots Mars's orbit in red.
        
        ax.plot(R(a_POD, e_POD, f_POD)*np.cos(f_POD), R(a_POD, e_POD, f_POD)*np.sin(f_POD), color='b', marker='o')              # Plots Earth's position upon TMI.
        ax.plot(R(a_target, e_target, f_target)*np.cos(f_target+Theta), R(a_target, e_target, f_target)*np.sin(f_target+Theta), color='r', marker='o')      # Plots Mars's position upon TMI.       
        
        ax.plot(R(a_POD, e_POD, f_PODOA)*np.cos(f_PODOA), R(a_POD, e_POD, f_PODOA)*np.sin(f_PODOA), color='b', marker='o', fillstyle='none')  # Plots Earth's position upon MOI.        
        ax.plot(R(a_target, e_target, f_targetOA)*np.cos(f_targetOA+Theta), R(a_target, e_target, f_targetOA)*np.sin(f_targetOA+Theta), color='g', marker='o', fillstyle='none')  # Plots Mars's position upon MOI.
        ax.plot(0, 0, color='k', marker='x') # Plots a black "x" to indicate the Sun's location.

        plt.savefig("Mission Profiles/result" + str(JDN) + ".eps", format="eps")        
        plt.show()
  
def plot_return(f_POD, f_PODOA, f_target, f_targetOA, Theta, a_trans, e_trans, a_POD, e_POD, a_target, e_target, JDN):
        fig = plt.figure() # initialize figure
        ax = fig.add_subplot(111) # name of the plot
        f = np.arange(0,2*np.pi,2*np.pi/999)
        ax.set_aspect('equal')
        
        R_T = a_trans*(1 - e_trans**2)/(1 + e_trans*np.cos(f))
        R_M = a_POD*(1 - e_POD**2)/(1 + e_POD*np.cos(f))
        R_E = a_target*(1 - e_target**2)/(1 + e_target*np.cos(f))
            
        ax.plot(R_T*np.cos(f+f_POD+Theta-np.pi), R_T*np.sin(f+f_POD+Theta-np.pi), color='y', lw=1)     # Plots the transfer's orbit in yellow.
        ax.plot(R_E*np.cos(f), R_E*np.sin(f), color='b', lw=1)             # Plots Earth's orbit in blue.
        ax.plot(R_M*np.cos(f+Theta), R_M*np.sin(f+Theta), color='r', lw=1) # Plots Mars's orbit in red.
        
        ax.plot(R(a_POD, e_POD, f_POD)*np.cos(f_POD+Theta), R(a_POD, e_POD, f_POD)*np.sin(f_POD+Theta), color='r', marker='o')              # Plots Mars's position upon TEI.
        ax.plot(R(a_target, e_target, f_target)*np.cos(f_target), R(a_target, e_target, f_target)*np.sin(f_target), color='b', marker='o')      # Plots Earth's position upon TEI.       
        
        ax.plot(R(a_POD, e_POD, f_PODOA)*np.cos(f_PODOA+Theta), R(a_POD, e_POD, f_PODOA)*np.sin(f_PODOA+Theta), color='r', marker='o', fillstyle='none')  # Plots Mars's position upon EOI.        
        ax.plot(R(a_target, e_target, f_targetOA)*np.cos(f_targetOA), R(a_target, e_target, f_targetOA)*np.sin(f_targetOA), color='y', marker='o', fillstyle='none')  # Plots Earth's position upon EOI.
        ax.plot(0, 0, color='k', marker='x') # Plots a black "x" to indicate the Sun's location.

        plt.savefig("Mission Profiles/result" + str(JDN) + ".eps", format="eps")        
        plt.show()
      
def plot_final(f_low, f_low_OA, f_low_return, f_high, f_high_OA, Theta, a_trans, e_trans, a_trans_return, e_trans_return, a_low, e_low, a_high, e_high, JDN):
        fig = plt.figure() # initialize figure
        ax = fig.add_subplot(111) # name of the plot
        f = np.arange(0,2*np.pi,2*np.pi/999)
        ax.set_aspect('equal')
        
        R_T = a_trans*(1 - e_trans**2)/(1 + e_trans*np.cos(f))
        R_T_return = a_trans_return*(1 - e_trans_return**2)/(1 + e_trans_return*np.cos(f))
        R_E = a_low*(1 - e_low**2)/(1 + e_low*np.cos(f))
        R_M = a_high*(1 - e_high**2)/(1 + e_high*np.cos(f))
            
        ax.plot(R_T*np.cos(f+f_low), R_T*np.sin(f+f_low), color='g', lw=1)     # Plots the transfer's orbit in green.
        ax.plot(R_T_return*np.cos(f+f_low), R_T_return*np.sin(f+f_low), color='y', lw=1)     # Plots the transfer's orbit in green.
        ax.plot(R_E*np.cos(f), R_E*np.sin(f), color='b', lw=1)             # Plots Earth's orbit in blue.
        ax.plot(R_M*np.cos(f+Theta), R_M*np.sin(f+Theta), color='r', lw=1) # Plots Mars's orbit in red.
        
        ax.plot(R(a_low, e_low, f_low)*np.cos(f_low), R(a_low, e_low, f_low)*np.sin(f_low), color='b', marker='o')              # Plots Earth's position upon TMI.
        #ax.plot(R(a_high, e_high, f_high)*np.cos(f_high+Theta), R(a_high, e_high, f_high)*np.sin(f_high+Theta), color='r', marker='o')      # Plots Mars's position upon TMI.       
        
        #ax.plot(R(a_low, e_low, f_low_OA)*np.cos(f_low_OA), R(a_low, e_low, f_low_OA)*np.sin(f_low_OA), color='b', marker='o', fillstyle='none')  # Plots Earth's position upon MOI.        
        ax.plot(R(a_high, e_high, f_high_OA)*np.cos(f_high_OA+Theta), R(a_high, e_high, f_high_OA)*np.sin(f_high_OA+Theta), color='g', marker='o', fillstyle='none')  # Plots Mars's position upon MOI.
        ax.plot(R(a_low, e_low, f_low_return)*np.cos(f_low_return), R(a_low, e_low, f_low_OA)*np.sin(f_low_OA), color='b', marker='o', fillstyle='none')  # Plots Earth's position upon MOI. 
        
        ax.plot(0, 0, color='k', marker='x') # Plots a black "x" to indicate the Sun's location.

        plt.savefig("result" + str(JDN) + ".eps", format="eps")        
        plt.show()

def plot_mission(mission_profiles, m, Theta):
    transfer_plot(mission_profiles[m,9], mission_profiles[m,13], mission_profiles[m,10], mission_profiles[m,12], Theta, mission_profiles[m,3], mission_profiles[m,4], mission_profiles[m,7], mission_profiles[m,8], mission_profiles[m,5], mission_profiles[m,6], mission_profiles[m,0])
    plot_return(mission_profiles[m,10+16], mission_profiles[m,13+16], mission_profiles[m,9+16], mission_profiles[m,12+16], Theta, mission_profiles[m,3+16], mission_profiles[m,4+16], mission_profiles[m,7+16], mission_profiles[m,8+16], mission_profiles[m,5+16], mission_profiles[m,6+16], mission_profiles[m,0+16])