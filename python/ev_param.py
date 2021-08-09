#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep  3 10:21:16 2020

@author: timowenzel
"""
import scipy.io 
import numpy as np
import random
import csv

def ev_load(path_input):

    # load ev data            
    ev_raw_bed = scipy.io.loadmat(path_input + "/ev_profiles_20cars_11kW_residential.mat")
    
    # convert 15 min data to hourly data
    ev = {}
    ev["avail"] = np.zeros(shape=(len(ev_raw_bed["EV_occurrence"][0]),8760))
    ev["dem_arrive"] = np.zeros(shape=(len(ev_raw_bed["EV_uncontrolled_charging_profile"][0]),8760))
    for n in range(len(ev_raw_bed["EV_occurrence"][0])):
        for t in range(8760):
            ev["avail"][n,t] = np.round(np.sum(ev_raw_bed["EV_occurrence"][t*4:(t*4+4),n])/4)  
    
    daily_dem = ev_raw_bed["EV_daily_demand"]      
    
    
    # map demand to end of available phase (grid-ractive and bi directional charging)
    index_leave = {}
    for i in range(len(ev["avail"])):
        index_leave[i] = ev["avail"][i,:-1] > ev["avail"][i,1:]
        # last charging phase of the year ends with t=T -> set last entry True
        index_leave[i] = np.append(index_leave[i], True)          
    

    # map demands to beginning of available phase (on-demand charging)
    index_arrive = {}
    index_temp = {}
    for i in range(len(ev["avail"])):
        index_temp[i] = ev["avail"][i,:-1] < ev["avail"][i,1:]
        index_arrive[i] = False
        index_arrive[i] = np.append(index_arrive[i], index_temp[i])
    
    
    ev["dem_leave"] = np.genfromtxt(open(path_input+"ev_dem_leave.csv", "rb"), delimiter = ",")

    
# =============================================================================
#     ev["dem_leave"] = np.zeros([8760,20])
#     for n in range(20):
#           for i in range(365):
#                kappa=random.uniform(7.5,15) 
#                ev["dem_leave"][(24*i):(24*i+24),n] = index_leave[n][(24*i):(24*i+24)] * kappa
#                ev["dem_leave"][8759,n] = index_leave[n][8759] * kappa 
#       
#       
#     np.savetxt("ev_dem_leave.csv", ev["dem_leave"], delimiter = ",")
# =============================================================================

    ev["dem_arrive"] = np.zeros([8760, 20])
    for n in range(20):
        for i in range(0,365):
            kappa = random.uniform(5,15)
            ev["dem_arrive"][(24*i):(24*i+24),n] = index_arrive[n][(24*i):(24*i+24)] * kappa
            
    return ev








