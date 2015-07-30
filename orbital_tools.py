# -*- coding: utf-8 -*-
"""
Created on Tue Jul 28 22:09:52 2015

@author: Mott
"""
from __future__ import division
import numpy as np


def mean_anomaly(e,M): 
    return

def eccentric_anomaly(e,M):    
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
    return (E1,count)

def eccentric_anomaly_old(e,M):    
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
    
def true_anomaly(e,E):   
    f = 2*np.arctan(((1 + e) / (1 - e))**(1/2)*np.tan(E/2)) # True Anomaly
    
    # Test to Make Sure True Anomaly is Between 0 and 360 Degrees
    if f < 0: 
        f = f + 2*np.pi
    elif f > 2*np.pi:
        f = f - 2*np.pi
    
    return f