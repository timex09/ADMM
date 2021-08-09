#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 11 16:20:57 2020
@author: timowenzel
"""
# import python packages
import numpy as np


# import own functions
import python.opti_stochastic_distributed as stoch_opti_dist

#%% GAME THEORY

def admm(nodes,devs,tech_par,settings):

    # Time Parameter
    time_hor = settings["time_horizon"]
    scenarios = settings["scenarios"]
    rounds = {}
    costs = {}
    res = {}
    optiparams = {}


    #%% Optimization Iteration
    status = False    # State of ADMM convergence
    counter = 0     # Iteration counter
    while status == False and counter < settings["admm_iteration_limit"]:

        if counter == 0:
            admm_params = { "lambda_load": settings["admm_startval_dual"] * np.ones_like(nodes["elec_scn"]),
                            "lambda_inj": settings["admm_startval_dual"] * np.ones_like(nodes["elec_scn"]),
                            "avg_load": nodes["diff_power_load"]/scenarios,
                            "avg_inj": nodes["diff_power_inj"]/scenarios
                            }

            for s in range(scenarios):
                costs[s], res[s], optiparams[s] = stoch_opti_dist.opti(nodes, devs, tech_par, settings, s, admm_params)

            # Safe results of each admm round
            rounds[counter] = res.copy()

        else:
            print("********************************")
            print("This is iteration No." + str(counter))
            print("********************************")

            for s in range(scenarios):
                costs[s], res[s], optiparams[s] = stoch_opti_dist.opti(nodes, devs, tech_par, settings, s, admm_params)

            # Safe results of each admm round
            rounds[counter] = res.copy()

            # Checking for convergence
            # todo: adjust to right temination criterion
            status_all = np.zeros(shape=[time_hor, scenarios])
            for t in range(time_hor):
                for s in range(scenarios):
                    if abs(rounds[counter - 1][s]["p_load_lem"][t] - rounds[counter][s]["p_load_lem"][t]) >= settings["addm_threshold"]:
                        status_all[t, s] = False
                    else:
                        status_all[t, s] = True

            status = status_all.all()

        # count up
        counter = counter + 1

        # param update
        avg_load = np.zeros_like(res[0]["p_load_lem"])
        avg_inj = np.zeros_like(res[0]["p_inj_lem"])
        for t in range(time_hor):
            avg_load[t] = sum(res[s]["p_load_lem"][t] for s in range(scenarios))/scenarios
            avg_inj[t] = sum(res[s]["p_inj_lem"][t] for s in range(scenarios))/scenarios

        lambda_load = np.zeros_like(nodes["elec_scn"])
        lambda_inj = np.zeros_like(nodes["elec_scn"])
        for t in range(time_hor):
            for s in range(scenarios):
                lambda_load[t, s] = admm_params["lambda_load"][t, s] + settings["admm_param_gamma"] * (
                        res[s]["p_load_lem"][t] - avg_load[t])
                lambda_inj[t, s] = admm_params["lambda_inj"][t, s] + settings["admm_param_gamma"] * (
                        res[s]["p_inj_lem"][t] - avg_load[t])


        admm_params = { "lambda_load": lambda_load,
                        "lambda_inj": lambda_inj,
                        "avg_load": avg_load,
                        "avg_inj": avg_inj
                        }

    return counter, status, costs, res, optiparams

