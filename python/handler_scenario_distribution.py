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

def admm(nodes, options, week, net_data, params):

    # Time Parameter
    time_hor = options["time_hor"]
    rounds = {}
    Speicher={}
    costs = {}
    res = {}
    optiparams = {}

    time_steps = list(range(0, options["time_hor"]))
    dt = 1
    last_time_step = options["time_hor"]
    idx = week


    #%% Optimization Iteration
    status = False    # State of ADMM convergence
    counter = 0     # Iteration counter

    AA=np.zeros((200, len(time_steps)))
    HilfsvariableA = AA
    HilfsvariableB = AA
    HilfsvariableC = AA
    HilfsvariableD = AA


    #for i in range(200):
        #HilfsvariableA[i] = 0
       # HilfsvariableB[i] = 0
        #HilfsvariableC[i] = 0
       # HilfsvariableD[i] = 0
        #for t in time_steps:
            #HilfsvariableA[i][t] = 0
            #HilfsvariableB[i][t] = 0
            #HilfsvariableC[i][t] =0
           # HilfsvariableD[i][t] = 0

    primal=np.zeros((len(nodes), len(time_steps)))
    dual=np.zeros((len(nodes), len(time_steps)))
    #for n in nodes:
        #primal[n]={}
        #dual[n]={}
       # for t in time_steps:
             #   primal[n][t]={}
               # dual[n][t] = {}

    #primal_Ges={}
    #dual_Ges={}
   # for i in range(200):
       # primal_Ges[i] = {}
        #dual_Ges[i]={}

    lam={}
    for n in nodes:
        lam[n]={}
        for t in time_steps:
            lam[n][t]={}


    #P_P2P_inj={}
    #P_P2P_load = {}
    #for n in nodes:
       # P_P2P_inj[n] = {}
       # P_P2P_load[n] = {}
       # for j in nodes:
           # P_P2P_inj[n][j] = {}
          #  P_P2P_load[n][j] = {}
           # for t in time_steps:
               # P_P2P_inj[n][j][t] = {}
               # P_P2P_load[n][j][t] = {}

    #BB = np.zeros((40,40, len(time_steps)+1))


    CC = np.zeros((200, len(nodes), len(nodes) , len(time_steps) ))
    Speicher["p_load_P2P"] = CC
    Speicher["p_inj_P2P"] = CC
    BB = np.zeros((len(nodes), len(nodes) , len(time_steps) ))
    P_P2P_inj=BB
    P_P2P_load=BB


    NN=np.zeros((200, len(nodes), len(time_steps)))
    HilfsvariableE = NN
    HilfsvariableF = NN

    Ubergabedaten = {"P_P2P_inj":BB,
                         "P_P2P_load":BB,
                         "P_Lade_Netz":HilfsvariableA,
                         "P_Einsp_Netz":HilfsvariableB,
                         "P_Nutz":HilfsvariableC,
                         "P_Verkauf":HilfsvariableD,
                        "P_Sammel_inj":HilfsvariableE,
                        "P_Sammel_load":HilfsvariableF}

    #np.zeros(shape=[time_hor, scenarios])

    while status == False and counter < options["admm_iteration_limit"]:
        print(counter)
        if counter == 0:
            print('yes')
            for n in nodes:
                for t in time_steps:
                    lam[n][t]=options["admm_startval_dual"]

            admm_params = { "lambda_load": lam
                            }


            for z in nodes:
                costs[z], res[z], optiparams[z] = stoch_opti_dist.opti(nodes, options, week, net_data, params, z, admm_params, Ubergabedaten, counter,lam)

            # Safe results of each admm round
                rounds[counter] = res.copy()
                for t in time_steps:
                    HilfsvariableA[counter][t] = HilfsvariableA[counter][t] + (res[z]["p_load_grid"][t])
                    HilfsvariableB[counter][t] = HilfsvariableB[counter][t] + (res[z]["p_inj_grid"][t])
                    HilfsvariableC[counter][t] = HilfsvariableC[counter][t] + (res[z]["p_use"][t])
                    HilfsvariableD[counter][t] = HilfsvariableD[counter][t] + (res[z]["p_sell"][t])
                    HilfsvariableE[counter][z][t] = res[z]["p_sell"][t]
                    HilfsvariableF[counter][z][t] = res[z]["p_use"][t]

                    for j in nodes:
                        P_P2P_inj[z][j][t] = res[z]["p_inj_P2P"][j][t]
                        P_P2P_load[z][j][t] = res[z]["p_load_P2P"][j][t]
                        Speicher["p_load_P2P"][counter][z][j][t] = P_P2P_load[z][j][t]
                        Speicher["p_inj_P2P"][counter][z][j][t] = P_P2P_inj[z][j][t]
                print(sum(res[z]["p_inj_P2P_Ges"][t] for t in time_steps))
                print(sum(res[z]["p_load_P2P_Ges"][t] for t in time_steps))
                print(sum(res[z]["p_load_grid"][t] for t in time_steps))

                Ubergabedaten = {"P_P2P_inj":P_P2P_inj,
                         "P_P2P_load":P_P2P_load,
                         "P_Lade_Netz":HilfsvariableA,
                         "P_Einsp_Netz":HilfsvariableB,
                         "P_Nutz":HilfsvariableC,
                         "P_Verkauf":HilfsvariableD,
                        "P_Sammel_inj": HilfsvariableE,
                        "P_Sammel_load": HilfsvariableF
                                 }
            for n in nodes:
                for t in time_steps:
                    primal[n][t] ==2
                    dual[n][t] ==2
                    print('angekommen')
            status=False
        else:
            print("********************************")
            print("This is iteration No." + str(counter))
            print("********************************")

            for z in nodes:
                costs[z], res[z], optiparams[z]= stoch_opti_dist.opti(nodes, options, week, net_data, params, z, admm_params, Ubergabedaten, counter,lam)
                print(counter)
            # Safe results of each admm round
                rounds[counter] = res.copy()
                for t in time_steps:
                    HilfsvariableA[counter][t] = HilfsvariableA[counter][t] + (res[z]["p_load_grid"][t] )
                    HilfsvariableB[counter][t] = HilfsvariableB[counter][t] + (res[z]["p_inj_grid"][t])
                    HilfsvariableC[counter][t] = HilfsvariableC[counter][t] + (res[z]["p_use"][t])
                    HilfsvariableD[counter][t] = HilfsvariableD[counter][t] + (res[z]["p_sell"][t])
                    HilfsvariableE[counter][z][t] = res[z]["p_sell"][t]
                    HilfsvariableF[counter][z][t] = res[z]["p_use"][t]

                    for j in nodes:
                        P_P2P_inj[z][j][t] = res[z]["p_inj_P2P"][j][t]
                        P_P2P_load[z][j][t] = res[z]["p_load_P2P"][j][t]
                        Speicher["p_load_P2P"][counter][z][j][t] = P_P2P_load[z][j][t]
                        Speicher["p_inj_P2P"][counter][z][j][t] = P_P2P_inj[z][j][t]



                Ubergabedaten = {"P_P2P_inj": P_P2P_inj,
                 "P_P2P_load": P_P2P_load,
                 "P_Lade_Netz": HilfsvariableA,
                 "P_Einsp_Netz": HilfsvariableB,
                 "P_Nutz": HilfsvariableC,
                 "P_Verkauf": HilfsvariableD,
                "P_Sammel_inj": HilfsvariableE,
                "P_Sammel_load": HilfsvariableF
                                 }
            for n in nodes:
                for t in time_steps:
                        primal[n][t]== abs(sum(P_P2P_load[n][j][t] - P_P2P_inj[j][n][t] for j in nodes))
                        dual[n][t]== abs(options["rho"] * sum((Speicher["p_load_P2P"][counter][n][j][t] - Speicher["p_load_P2P"][counter-1][n][j][t] + Speicher["p_inj_P2P"][counter][n][j][t] - Speicher["p_inj_P2P"][counter-1][n][j][t]) for j in nodes))
            print(sum(primal[1][t] for t in time_steps))

            status_all = np.full((len(nodes), len(time_steps)), False, dtype=bool)
            for n in nodes:
                for t in time_steps:
                    if primal[n][t] > options["admm_threshold"] and dual[n][t] > options["admm_threshold"]:
                        # abs(rounds[counter - 1][s]["p_load_lem"][t] - rounds[counter][s]["p_load_lem"][t]) >= settings["addm_threshold"]:
                        status_all[n, t] = False
                    else:
                        status_all[n, t] = True

            status = status_all.all()

            # Checking for convergence
            # todo: adjust to right temination criterion
        #status_all = np.zeros(shape=[time_hor])
        #if counter == 0:
            #for n in nodes:
               # for t in time_steps:
                 #   primal[n][t] ==2
                  #  dual[n][t] ==2
                   # print('angekommen')
        #else:
          #  for n in nodes:
               # for t in time_steps:
                   #     primal[n][t]== abs(sum(P_P2P_load[n][j][t] - P_P2P_inj[j][n][t] for j in nodes))
                     #   dual[n][t]== abs(options["rho"] * sum((Speicher["p_load_P2P"][counter][n][j][t] - Speicher["p_load_P2P"][counter-1][n][j][t] + Speicher["p_inj_P2P"][counter][n][j][t] - Speicher["p_inj_P2P"][counter-1][n][j][t]) for j in nodes))

        #primal_Ges[counter]=sum(primal[n][t] for n in nodes for t in time_steps)
        #dual_Ges[counter] = sum(dual[n][t] for n in nodes for t in time_steps)
        #status_all = np.full((len(nodes), len(time_steps)), False, dtype=bool)

        #for n in nodes:
           # for t in time_steps:
                #if primal[n][t]>options["admm_threshold"] and dual[n][t]>options["admm_threshold"]:
                            #abs(rounds[counter - 1][s]["p_load_lem"][t] - rounds[counter][s]["p_load_lem"][t]) >= settings["addm_threshold"]:
                    #status_all[n, t] = False
               # else:
                   #status_all[n, t] = True

        #status = status_all.all()
        #
        # count up

        counter = counter + 1

        # param update
        print(status)

        #lambda_load = np.zeros_like(nodes)
        #lambda_inj = np.zeros_like(nodes)
        for t in time_steps:
            for n in nodes:
                lam[n][t] = lam[n][t] + options["rho"] * (primal[n][t])
                #lambda_inj[t, s] = admm_params["lambda_inj"][t, s] + settings["admm_param_gamma"] * (
                     #   res[s]["p_inj_lem"][t] - avg_load[t])




        admm_params = { "lambda_load": lam

                        }
        print(counter)

    print('ENDE')
    result_file = {}
    result_file["1"] = [lam[1][t] for t in time_steps]
    result_file["2"] = [lam[2][t] for t in time_steps]
    result_file["7"] = [lam[7][t] for t in time_steps]
    return counter, status, costs, res, optiparams, rounds, result_file



