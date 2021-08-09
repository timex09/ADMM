#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 03 2021

@author: she
"""
# import python packages
import numpy as np

# import own functions
#import python.parameters as parameters
#import python.building_handler as handler

import python.handler_scenario_distribution as verteilt
import python.load_params as parameters
import python.load_demands as demand
import python.load_devices as devices
import pandas as pd


def get_opti_inputs(options):

    # Load Params
    devs_temp = devices.load_devices()
    nodes_temp, building_params = demand.load_demands(options, devs_temp)
    nodes = devices.map_devices(nodes_temp, building_params, devs_temp)
    params = parameters.load_eco()
    net_data = parameters.load_net()


    return nodes, net_data, params


# %% GENERAL PARAMETERS
if __name__ == "__main__":



    options = {"time_hor": 168,
               "tweeks": 4,
               "buildings": 20,
               "pv_share": 0.5,
               "ev_status": True,
               "district": True,  # True -> district optimisation, False -> building optimization
               "grid": False,  # True -> consider grid constraints, False -> dont
               "rho": 1,
               "admm_startval_dual": (parameters.load_eco()["eco"]["elec_price"] + parameters.load_eco()["eco"]["eeg_pv"]) / 2,
               "admm_iteration_limit": 200,
               "admm_threshold": 0.1,
               "path_input": "C:/Users/Tim/Desktop/Ma EBC/Python/ADMM/input_data/bedburg/"}



    nodes, net_data, params = get_opti_inputs(options)

    for week in range(options["tweeks"]):
        counter, status, costs, res, optiparams, rounds, result = verteilt.admm(nodes, options, week, net_data, params)


    #with pd.ExcelWriter('output.xlsx') as writer:

           # result["from_grid"].to_excel(writer, sheet_name=str(i))


    print(result["1"])
    print(result["2"])
    print(result["7"])





    # %% Store the results into a pkl-date
    #import pickle

    #filename = "results/results_" + "mode" + str(settings["mode"]) + ".pkl"
   # with open(filename, "wb") as f_in:
       # pickle.dump(costs, f_in, pickle.HIGHEST_PROTOCOL)
        #pickle.dump(res, f_in, pickle.HIGHEST_PROTOCOL)
        #pickle.dump(nodes, f_in, pickle.HIGHEST_PROTOCOL)
        #pickle.dump(optiparams, f_in, pickle.HIGHEST_PROTOCOL)
       # pickle.dump(settings, f_in, pickle.HIGHEST_PROTOCOL)
       # pickle.dump(tech_par, f_in, pickle.HIGHEST_PROTOCOL)
