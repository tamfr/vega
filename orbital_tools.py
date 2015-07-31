# -*- coding: utf-8 -*-
"""
Created on Tue Jul 28 22:09:52 2015

@author: Mott
"""
from __future__ import division
import numpy as np


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
    
#######################

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
    
def t_of_M(a,mu,M):
    """Mean anomaly as fucntion of eccentric anomaly.
        Requires semi-major axis and standard gravitational parameter of the central body.
    """ 
    return M*(a**3/mu)**(1/2)