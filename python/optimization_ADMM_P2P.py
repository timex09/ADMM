#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Juni  24 07:48:52 2021

@author: Tim Adamek
"""

import numpy as np
import gurobipy as gp
import time
import pandas as pd
import python.store_items as store

def compute(nodes, options, week, net_data, params):

    now = time.time()

    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    # Setting up the model

    # Create a new model
    model = gp.Model("Operational_Optimization")

    time_steps = list(range(0, options["time_hor"]))
    dt = 1
    last_time_step = options["time_hor"]
    idx = week
    rho=1  #ADMM Parameter
    it=list(range(1, 200))
    eps_prim=0.01
    eps_dual=0.01
    lam = {}
    for n in nodes:
        lam[n] = {}
        for t in time_steps:
            lam[n][t] = {}
            for p in  it:
                lam[n][t][p]={}

# Startwert Lambda
    for n in nodes:
        for t in time_steps:
            lam[n][t][1]=(params["eco"]["elec_price"]+params["eco"]["eeg_pv"])/2

    for t in time_steps:
        for i in it:
            while primal_Ges[t][i]>eps_prim & dual_Ges[t][i]>eps_dual:


# Zielfunktion je Haushalt
    model.addConstr(
        operational_costs == power_grid["load"][n][t]*params["eco"]["elec_price"] - power_grid["inj"][n][t]*params["eco"]["eeg_pv"] + (sum(power_P2P["load"][n][j][t]*lam[j][t][??]) for j in nodes)
        + power_P2P_sum["inj"][n][t]*lam[n][t][??] + rho/2 * (sum(((power_P2P["inj"][n][j][t]-power_P2P["load"][j][n][t])**2)+((power_P2P["inj"][j][n][t]-power_P2P["load"][n][j][t])**2)) for j in nodes)

# Residuen
    primal = {}
    for n in nodes:
        primal[n] = {}
        for t in time_steps:
            primal[n][t] = {}
            for p in it:
                primal[n][t][p] = {}

    primal[n][t][p] == sum((power_P2P["inj"][n][j][t]-power_P2P["load"][j][n][t]) for j in nodes)

    primal_Ges=[]
    for t in time_steps:
        primal_Ges[t]= {}
        for p in it_
            primal_Ges[t][p]=sum(primal[n][t][p] for n in nodes)

    dual = {}
    for n in nodes:
        dual[n] = {}
        for t in time_steps:
            dual[n][t] = {}
            for p in it:
                dual[n][t][p] = {}
    dual[n][t][p] == rho * (power_P2P["inj"][n][j][t][p]-power_P2P["inj"][j][n][t][p-1] + power_P2P["load"][n][j][t][p]-power_P2P["load"][j][n][t][p-1])

    dual_Ges = []
    for t in time_steps:
        dual_Ges[t] = {}
        for p in it_
            dual_Ges[t][p]=sum(dual[n][t][p] for n in nodes)

    # Lambda Update
    lam[n][t][p+1] == lam[n][t][p] + rho * primal[n][t][p]

    # %% TECHNICAL PARAMETERS
    soc_nom = {}
    soc_nom["TES"] = {}
    soc_nom["BAT"] ={}
    soc_nom["EV"] = {}
    soc_init = {}
    soc_init["TES"] = {}
    soc_init["BAT"] = {}
    soc_init["EV"] = {}
    for n in nodes:
        # Define nominal SOC_nom according to capacities: SOC: Ladezustand 0,1=10%
        soc_nom["TES"][n] = params["phy"]["c_w"] * nodes[n]["TES"]["cap"] * \
                    nodes[n]["TES"]["dT_max"] * params["phy"]["rho_w"] / 3600
        soc_nom["BAT"][n] = nodes[n]["BAT"]["cap"]
        if options["ev_status"]:
            soc_nom["EV"][n]= nodes[n]["EV"]["cap"]  # kWh   Nomincal capcity Battery from ev
        else:
            soc_nom["EV"][n] = 0
        # Initial state of charge
        soc_init["TES"][n] = soc_nom["TES"][n] * 0.5  # kWh   Initial SOC TES
        soc_init["BAT"][n] = soc_nom["BAT"][n] * 0.5  # kWh   Initial SOC Battery
        soc_init["EV"][n] = soc_nom["EV"][n] * 0.75

        # %% TECHNICAL CONSTRAINTS

    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    # CREATE VARIABLES

    #%% OPERATIONAL BUILDING VARIABLES
    #   NOTE: Subscript "dom" = "domestic"
    # Electrical power to/from electricity-based domestic devices
    power_dom = {}
    for device in ["HP", "EB", "PV"]:
        power_dom[device] = {}
        for n in nodes:
            power_dom[device][n] = {}
            for t in time_steps:
                power_dom[device][n][t] = model.addVar(vtype="C", name="power_" + device + "_n" + str(n) + "_t" + str(t))


    # Heat to/from devices
    heat_dom = {}
    for device in ["HP", "EB"]:
       heat_dom[device] = {}
       for n in nodes:
           heat_dom[device][n] = {}
           for t in time_steps:
               heat_dom[device][n][t] = model.addVar(vtype="C", name="heat_"  + device + "_n" + str(n) + "_t" + str(t))

    # Storage variables
    soc_dom = {}                # State of charge
    ch_dom = {}
    dch_dom = {}
    for device in ["TES","BAT","EV"]:
        soc_dom[device] = {}    #Energy (kWh)
        ch_dom[device] = {}     #Power(kW)
        dch_dom[device] = {}    #Power (kW)
        for n in nodes:
            soc_dom[device][n] = {}
            ch_dom[device][n] = {}
            dch_dom[device][n] = {}
            for t in time_steps:
                soc_dom[device][n][t] = model.addVar(vtype="C", name="soc_" + device + "_n" + str(n) + "_t" + str(t))
                ch_dom[device][n][t] = model.addVar(vtype="C", name="ch_dom_" + device + "_n" + str(n) + "_t" + str(t))
                dch_dom[device][n][t] = model.addVar(vtype="C", name="dch_dom_" + device + "_n" + str(n) + "_t" + str(t))

    # Residual building demands (in kW)  [Sum for each building of all devices]
    res_dom = {}
    res_dom["power"] = {}
    res_dom["feed"] = {}
    for n in nodes:
        # Electricity demand
        res_dom["power"][n] = {}
        res_dom["feed"][n] = {}
        for t in time_steps:
            res_dom["power"][n][t] = model.addVar(vtype="C", name="residual_power_n" + str(n) + "_t" + str(t))
            res_dom["feed"][n][t] = model.addVar(vtype="C", name="residual_feed_n" + str(n) + "_t" + str(t))

    # P2P Handelsmengen: Netz
    power_grid={}
    power_grid["inj"] ={}
    power_grid["load"] = {}
    for n in nodes:
       power_grid["inj"][n] = {}
       power_grid["load"][n] = {}
       for t in time_steps:
            power_grid["inj"][n][t] = model.addVar(vtype="C", name="power_from_grid" + str(n) + "_t" + str(t))
            power_grid["load"][n][t] = model.addVar(vtype="C", name="power_to_grid" + str(n) + "_t" + str(t))


    # P2P Handelsmengen: P2P Leistung
    power_P2P={}
    power_P2P["inj"] = {}
    power_P2P["load"] = {}
    for n in nodes:
        power_P2P["inj"][n] = {}
        power_P2P["load"][n] = {}
        for j in nodes:
            power_P2P["inj"][n][j] = {}
            power_P2P["load"][n][j] = {}
            for t in time_steps:
                power_P2P["inj"][n][j][t] = model.addVar(vtype="C", name="power_from_P2P" + str(n) + "_t" + str(t))
                power_P2P["load"][n][j][t] = model.addVar(vtype="C", name="power_to_P2P" + str(n) + "_t" + str(t))



    # P2P Handelsmengen: P2P Leistung
    power_P2P_sum={}
    power_P2P_sum["inj"] = {}
    power_P2P_sum["load"] = {}
    for n in nodes:
        power_P2P_sum["inj"][n] = {}
        power_P2P_sum["load"][n] = {}
        for t in time_steps:
                power_P2P_sum["inj"][n][t] = model.addVar(vtype="C", name="power_from_P2P_sum" + str(n) + "_t" + str(t))
                power_P2P_sum["load"][n][t] = model.addVar(vtype="C", name="power_to_P2P_sum" + str(n) + "_t" + str(t))




    #binary variable for each house to avoid simuultaneous feed-in and purchase of electric energy
    #binary = {}
   # for device in ["HLINE","BAT","EV"]:
       # binary[device] = {}
       # for n in nodes:
          #  binary[device][n] = {}
           # for t in time_steps:
               # binary[device][n][t] = model.addVar(vtype= gp.GRB.BINARY, name="factor_binary_" + device + "_n" + str(n) + "_t" + str(t))

    # Residual network demand
    residual = {}
    residual["power"] = {}          # Residual network electricity demand
    residual["feed"] = {}           # Residual feed in
    for t in time_steps:
        residual["power"][t] = model.addVar(vtype = "C", name="residual_power_t" + str(t))
        residual["feed"][t] = model.addVar(vtype = "C", name="residual_feed_t" + str(t))


    #%% NETWORK VARIABLES

    # Network variables
    power_net = {}
    power_net["Inj"] = {}
    power_net["Load"] = {}
    for node in net_data["gridnodes"]:
        # Electricity demand
        power_net["Inj"][node] = {}
        power_net["Load"][node] = {}
        for t in time_steps:
            power_net["Inj"][node][t] = model.addVar(vtype="C", name="powerInj_" + str(n) + "_t" + str(t))
            power_net["Load"][node][t] = model.addVar(vtype="C", name="powerLoad_" + str(n) + "_t" + str(t))

    # set line bounds due to technical limits
    powerLine = model.addVars(net_data["nodeLines"],time_steps, vtype="C", lb=-10000, name="powerLine_")

    # set trafo bounds due to technichal limits
    powerTrafoLoad = model.addVars(time_steps, vtype="C", lb=0, ub=net_data["trafo_max"], name="powerTrafo_"+str(t))
    powerTrafoInj = model.addVars(time_steps, vtype="C", lb=0, ub=net_data["trafo_max"], name="injTrafo_"+str(t))

    # activation variable for trafo load
    #yTrafo = model.addVars(time_steps, vtype="B", name="yTrafo_"+str(t))


    #%% BALANCING UNIT VARIABLES

    # Electrical power to/from devices
    power = {}
    for device in ["from_grid", "to_grid", "PV"]:
        power[device] = {}
        for t in time_steps:
            power[device][t] = model.addVar(vtype="C", lb = 0, name="power_" + device + "_t" + str(t))

    # total energy amounts taken from grid
    from_grid_total_el = model.addVar(vtype = "C", name= "from_grid_total_el")
    # total power to grid
    to_grid_total_el = model.addVar(vtype = "C", name="to_grid_total_el")

    # Total gross CO2 emissions
    co2_total = model.addVar(vtype = "c", lb=-gp.GRB.INFINITY, name = "total_CO2")

    # Peak network load
    network_load = model.addVar(vtype = "c", lb=-gp.GRB.INFINITY, name = "peak_network_load")

    # Total operational costs
    operational_costs = model.addVar(vtype = "C", lb=-gp.GRB.INFINITY, name = "total_operational_costs")

    # Objective function
    obj = model.addVar(vtype="C", lb=-gp.GRB.INFINITY, name="obj")


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    # DEFINE OBJECTIVE FUNCTION
    model.update()
    model.setObjective(obj)
    model.ModelSense = gp.GRB.MINIMIZE


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    # ADD CONSTRAINTS

    # Device generation <= device capacity
    for n in nodes:
        for t in time_steps:
            for device in ["HP", "EB"]:
                model.addConstr(heat_dom[device][n][t] <= nodes[n][device]["cap"], name = str(device) + "_cap_" + str(n))

    for n in nodes:
        for t in time_steps:
            # Electric boiler
            model.addConstr(heat_dom["EB"][n][t] == power_dom["EB"][n][t], name = "EB_energybalance_" + str(n) + "_" + str(t))
            # Heat Pump
            model.addConstr(heat_dom["HP"][n][t] == power_dom["HP"][n][t] * nodes[n]["COP_sh"][:,idx][t], name = "heatpump_energybalance_heating" + str(n) + "_" + str(t))

    for n in nodes:
        for t in time_steps:
            model.addConstr(power_dom["PV"][n][t] <= nodes[n]["pv_power"][:,idx][t], name = "max_pv_power_" + str(n) + "_" + str(t))


    # SOC of storages <= storage capacity
    for n in nodes:
        for device in ["TES"]:
            for t in time_steps:
                model.addConstr(soc_dom[device][n][t] <= soc_nom["TES"][n], name ="max_cap_" +str(device)+ "_" + str(n) + "_" + str(t))
                model.addConstr(soc_dom[device][n][t] >= nodes[n]["TES"]["min_soc"] * soc_nom["TES"][n], name ="min_cap_" + str(device)+ "_" + str(n) + "_" + str(t))

    for n in nodes:
        for t in time_steps:
            model.addConstr(ch_dom["BAT"][n][t] <=  nodes[n]["BAT"]["max_ch"], name ="max_ch_cap_bat_" + str(n) + "_" + str(t))
            model.addConstr(dch_dom["BAT"][n][t] <=  nodes[n]["BAT"]["max_dch"], name ="max_dch_cap_bat_" + str(n) + "_" + str(t))

            #model.addConstr(ch_dom["BAT"][n][t] <= binary["BAT"][n][t]*100, name ="Binary1__bat" + str(n) + "_" + str(t))
            #model.addConstr(dch_dom["BAT"][n][t] <= (1-binary["BAT"][n][t])*100, name ="Binary2__bat" + str(n) + "_" + str(t))


            model.addConstr(soc_dom["BAT"][n][t] <=  nodes[n]["BAT"]["max_soc"] * nodes[n]["BAT"]["cap"], name ="max_soc_bat_" + str(n) + "_" + str(t))
            model.addConstr(soc_dom["BAT"][n][t] >=  nodes[n]["BAT"]["min_soc"] * nodes[n]["BAT"]["cap"], name ="max_soc_bat_" + str(n) + "_" + str(t))

#P2P Constraints
    for n in nodes:
        for j in nodes:
            for t in time_steps:
                if n==j:
                    model.addConstr(power_P2P["inj"][n][j][t]==0, name ="P2P_Constraint" + str(n) + "_" + str(j)+ "_" + str(t))
                    model.addConstr(power_P2P["load"][n][j][t] == 0, name ="P2P_Constraint" + str(n) + "_" + str(j)+ "_" + str(t))
                #else:
                    #model.addConstr(power_P2P["inj"][n][j][t]==power_P2P["load"][j][n][t], name ="P2P_Constraint" + str(n) + "_" + str(j)+ "_" + str(t))
                    #model.addConstr(power_P2P["inj"][j][n][t] == power_P2P["load"][n][j][t])


    for n in nodes:
        for t in time_steps:

            model.addConstr(power_P2P_sum["inj"][n][t] == sum(power_P2P["inj"][n][j][t] for j in nodes), name ="P2P" + str(n) + "_" + str(j)+ "_" + str(t))
            model.addConstr(power_P2P_sum["load"][n][t] == sum(power_P2P["load"][n][j][t] for j in nodes),name ="P2P" + str(n) + "_" + str(j)+ "_" + str(t))

    for n in nodes:
            for t in time_steps:
                model.addConstr(res_dom["power"][n][t]==power_P2P_sum["load"][n][t] + power_grid["load"][n][t])
                model.addConstr(res_dom["feed"][n][t]==power_P2P_sum["inj"][n][t] + power_grid["inj"][n][t])

    for t in time_steps:
        model.addConstr(sum(power_grid["load"][n][t] for n  in nodes) ==power["from_grid"][t])
        model.addConstr(sum(power_grid["inj"][n][t] for n in nodes) == power["to_grid"][t])

   # for t in time_steps:
       # model.addConstr((sum(power_P2P_sum["load"][n][t] - power_P2P_sum["inj"][n][t])for n in nodes)==0)





#%% EV CONSTRAINTS CONSTRAINTS
    device = "EV"
    if options["ev_status"]:
        for n in nodes:
            # Energy balance
            for t in time_steps:
                if t == 0:
                    soc_prev = soc_init[device][n]

                else:
                    soc_prev = soc_dom[device][n][t-1]

                model.addConstr(soc_dom[device][n][t] == soc_prev
                                        + ch_dom[device][n][t] * nodes[n][device]["eta_ch_ev"] * dt
                                        - dch_dom[device][n][t] / nodes[n][device]["eta_dch_ev"] * dt - nodes[n]["ev_dem_leave"][t])


                if t == last_time_step:

                    model.addConstr(soc_dom[device][n][t] == soc_init[device][n])


                model.addConstr(dch_dom[device][n][t] <= nodes[n]["ev_avail"][t] * nodes[n][device]["max_dch_ev"])
                model.addConstr(ch_dom[device][n][t]  <= nodes[n]["ev_avail"][t] * nodes[n][device]["max_ch_ev"])

                #model.addConstr(dch_dom[device][n][t] <= binary[device][n][t]*nodes[n][device]["max_dch_ev"], name ="Binary1_ev" + str(n) + "_" + str(t))
                #model.addConstr(ch_dom[device][n][t] <= (1-binary[device][n][t])*nodes[n][device]["max_ch_ev"], name ="Binary2_ev" + str(n) + "_" + str(t))

                model.addConstr(soc_dom[device][n][t] >= nodes[n][device]["min_soc"] * nodes[n][device]["cap"])
                model.addConstr(soc_dom[device][n][t] <= nodes[n][device]["max_soc"] * nodes[n][device]["cap"])

    else:
        for n in nodes:
            for t in time_steps:
                model.addConstr(soc_dom[device][n][t] == 0)
                model.addConstr(dch_dom[device][n][t] == 0)
                model.addConstr(ch_dom[device][n][t] == 0)


    #%% DOMESTIC FLEXIBILITIES

    # SOC coupled over all times steps (Energy amount balance, kWh)
    for n in nodes:
        # Energy balance TES
        device = "TES"
        for t in time_steps:
            if t == 0:
                soc_prev = soc_init[device][n]
            else:
                soc_prev = soc_dom[device][n][t - 1]

            model.addConstr(soc_dom[device][n][t] == soc_prev * nodes[n][device]["eta_tes"]
                                + ch_dom[device][n][t] * dt
                                - dch_dom[device][n][t] * dt, name ="TES_storage_balance_" + str(n) + "_" + str(t))

            model.addConstr(ch_dom[device][n][t] == heat_dom["HP"][n][t], name="Heat_charging_"+str(n) + "_" + str(t))
            model.addConstr(dch_dom[device][n][t] <= nodes[n]["heat"][:,idx][t] + 0.5 * nodes[n]["dhw"][:,idx][t], name="Heat_discharging_"+str(n) + "_" + str(t))

            if t == last_time_step:
                    model.addConstr(soc_dom[device][n][t] == soc_init[device][n], name="End_TES_Storage_"+str(n) + "_" + str(t))

        # Energy balance battery
        device = "BAT"
        for t in time_steps:
            if t == 0:
                soc_prev = soc_init[device][n]
            else:
                soc_prev = soc_dom[device][n][t - 1]

            model.addConstr(soc_dom[device][n][t] == soc_prev
                                    + ch_dom[device][n][t] * nodes[n][device]["eta_bat"] * dt
                                    - dch_dom[device][n][t] / nodes[n][device]["eta_bat"] * dt, name ="BAT_storage_balance_" + str(n) + "_" + str(t))

            if t == last_time_step:
                    model.addConstr(soc_dom[device][n][t] == soc_init[device][n], name="End_BAT_Storage_"+str(n) + "_" + str(t))

    #%% ENERGY BALANCES (Power balance, kW)

    for n in nodes:
        for t in time_steps:
            model.addConstr(res_dom["power"][n][t] + power_dom["PV"][n][t] + dch_dom["BAT"][n][t] + dch_dom["EV"][n][t]
                        == nodes[n]["elec"][:,idx][t]  +power_dom["HP"][n][t] + power_dom["EB"][n][t] + ch_dom["BAT"][n][t] + ch_dom["EV"][n][t] + res_dom["feed"][n][t], name ="Electricity_balance_" + str(n) + "_" + str(t))

            model.addConstr(res_dom["feed"][n][t] <= power_dom["PV"][n][t] + dch_dom["BAT"][n][t] + dch_dom["EV"][n][t], name ="Feed-in_max_" + str(n) + "_" + str(t))


    for n in nodes:
        for t in time_steps:
            # Heat balance
            model.addConstr(heat_dom["EB"][n][t] == 0.5 * nodes[n]["dhw"][:,idx][t], name ="Heat_balance_" + str(n) + "_" + str(t))


  #  for n in nodes:
       # for t in time_steps:
           # model.addConstr(res_dom["power"][n][t]<= binary["HLINE"][n][t]*100, name ="Binary1_" + str(n) + "_" + str(t))
           # model.addConstr(res_dom["feed"][n][t]<= (1-binary["HLINE"][n][t])*100, name ="Binary2_" + str(n) + "_" + str(t))

    #%% SUM UP: BUILDING CONSTRAINTS

    # Residual loads
    for t in time_steps:
        # Residual network electricity demand (Power balance, MW)
        model.addConstr(residual["power"][t] == sum(res_dom["power"][n][t] for n in nodes))
        model.addConstr(residual["feed"][t] == sum(res_dom["feed"][n][t] for n in nodes))

    #%% ENERGY BALANCES: NETWORK AND ENERGY HUB

    for t in time_steps:

        # For all modes and scenarios:
        # Electricity balance (Power balance, MW)
        model.addConstr(residual["feed"][t] + power["from_grid"][t] == residual["power"][t] + power["to_grid"][t])

        model.addConstr(power["to_grid"][t] <= residual["feed"][t])



        #model.addConstr(power["from_grid"][t] <= yTrafo[t]*1000, name ="Binary1_" + str(n) + "_" + str(t))
        #model.addConstr(power["to_grid"][t]<= (1-yTrafo[t])*1000, name ="Binary2_" + str(n) + "_" + str(t))

    #%% NETWORK CONSTRAINTS

    if options["grid"]:

        for t in time_steps:

            power["from_grid"][t] == powerTrafoLoad[t]
            power["to_grid"][t] == powerTrafoInj[t]

        for t in time_steps:

            for index, row in net_data["grid_allo"].iterrows():

                model.addConstr(power_net["Inj"][row["gridnodes"]][t] == res_dom["feed"][row["nodes"]][t])
                model.addConstr(power_net["Load"][row["gridnodes"]][t] == res_dom["power"][row["nodes"]][t])

        for node in net_data["gridnodes"]:

            for t in time_steps:

                if node in net_data["net_nodes"]["load"] :

                    dummy = 1

                else:

                    model.addConstr(power_net["Inj"][node][t] == 0)
                    model.addConstr(power_net["Load"][node][t] == 0)

        for node in net_data["gridnodes"]:

            for t in time_steps:

                if node in net_data["net_nodes"]["trafo"]:

                    model.addConstr(powerLine.sum(node,'*',t) - powerLine.sum('*',node,t) ==
                        powerTrafoLoad[t] - powerTrafoInj[t], name="node balance_"+str(n))

                    model.addConstr(power_net["Inj"][node][t] == 0)
                    model.addConstr(power_net["Load"][node][t] == 0)

                else:

                    model.addConstr(powerLine.sum(node,'*',t) - powerLine.sum('*',node,t) ==
                                    power_net["Inj"][node][t] - power_net["Load"][node][t], name="node balance_"+str(node))

        for [n,m] in net_data["nodeLines"]:

                    for t in time_steps:
                        model.addConstr(powerLine[n,m,t] <= net_data["powerLine_max"][n,m], name="line power max_"+str(n)+str(m)+str(t))
                        model.addConstr(powerLine[n,m,t] >= (-1)*net_data["powerLine_max"][n,m], name="line power min_"+str(n)+str(m)+str(t))

    else:
        pass

    #%% SUM UP TOTAL ELECTRICITY FROM AND TO GRID

    # Total electricity amounts taken from grid (Energy amounts, MWh)
    model.addConstr(from_grid_total_el == sum(power["from_grid"][t] for t in time_steps))

    # Total electricity feed-in (Energy amounts, MWh)
    model.addConstr(to_grid_total_el == sum(power["to_grid"][t] for t in time_steps))





    #%% OBJECTIVE FUNCTIONS
    #select the objective function based on input parameters

    if options["district"]:
        # Total operational costs
        model.addConstr(
            operational_costs == from_grid_total_el * params["eco"]["elec_price"]
            + network_load * params["eco"]["peak_price"] - to_grid_total_el * params["eco"]["eeg_pv"])
        # Total operational costs without peak price
        #model.addConstr(
        #    operational_costs == from_grid_total_el * params["eco"]["elec_price"]
         #   - to_grid_total_el * params["eco"]["eeg_pv"])

        # Emissions
        model.addConstr(
            co2_total == from_grid_total_el / 1000 * params["eco"]["co2_el"])
        # Network load
        model.addConstrs(
            network_load >= power["from_grid"][t] for t in time_steps)


    else:
        # Users operational costs
        model.addConstr(
            operational_costs == sum(
                (power_grid["load"][n][t] * params["eco"]["elec_price"] - power_grid["inj"][n][t] * params["eco"]["eeg_pv"]- power_P2P_sum["inj"][n][t] * c_P2P[n][t] + power_P2P["load"][n][j][t] * c_P2P[j][t]) for t in time_steps for n in nodes for j in nodes))

        # Emissions
        model.addConstr(
            co2_total == sum((res_dom["power"][n][t] / 1000 * params["eco"]["co2_el"]) for t in time_steps for n in nodes))

    # objective
    model.addConstr(obj == operational_costs)
    #model.params.NonConvex = 2
    # Carry out optimization
    for t in time_steps:
        for i in it:
            while primal_Ges[t][i] > eps_prim & dual_Ges[t][i] > eps_dual:
                for n in nodes:
                    model.optimize()

                    primal[n][t][i] == sum((power_P2P["inj"][n][j][t] - power_P2P["load"][j][n][t]) for j in nodes)
                    dual[n][t][i] == rho * (power_P2P["inj"][n][j][t][i] - power_P2P["inj"][j][n][t][i - 1] +
                                            power_P2P["load"][n][j][t][i] - power_P2P["load"][j][n][t][i - 1])
                    # Lambda Update
                    lam[n][t][i + 1] == lam[n][t][i] + rho * primal[n][t][i]
            dual_Ges[t][i] = sum(dual[n][t][i] for n in nodes)
            primal_Ges[t][p] = sum(primal[n][t][p] for n in nodes)
            else: break
   # model.write('test.lp')


    later = time.time()
    difference = later - now
    print("********************************************")
    print("Model run time was " + str(difference) + " seconds")
    print("********************************************")

    if model.status==gp.GRB.Status.INFEASIBLE or model.status==gp.GRB.Status.INF_OR_UNBD:
        model.computeIIS()
        f=open('errorfile_hp.txt','w')
        f.write('\nThe following constraint(s) cannot be satisfied:\n')
        for c in model.getConstrs():
            if c.IISConstr:
                f.write('%s' % c.constrName)
                f.write('\n')
        f.close()


    #%% SAVE RESULTS IN ONE CENTRAL RESULT FILE: result_file
    r={}
    r["pi"] = {}
    for n in nodes:
        r["pi"][n] = {}
        for j in nodes:
            r["pi"][n][j] = {}
            for t in time_steps:
                r["pi"][n][j][t] ={}


    for n in nodes:
        for j in nodes:
            for t in time_steps:
                r["pi"][n][j][t]=(model.getConstrByName("P2P_Constraint" + str(n) + "_" + str(j)+ "_" + str(t))).getAttr('Pi')

    #for c in model.getConstrs():
       # print("Name" + str(c.constrName) + " Dual"+ str(c.pi))


    result_file = {}
    result_file["CO2_total"] = co2_total.X
    result_file["costs_total"] = operational_costs.X
    result_file["peak"] = network_load.X
    result_file["from_grid_total_el"] = from_grid_total_el.X
    result_file["to_grid_total_el"] = to_grid_total_el.X

    #residual power
    result_file["residual_power"] = pd.DataFrame()
    for i in range(len(nodes)):
        temp = [res_dom["power"][i][t].X for t in time_steps]
        result_file["residual_power"].insert(i, "Building_" + str(i), temp, True)

    #power_P2P_sum["inj"][n][t]
    test=[r["pi"][2][3][t] for t in time_steps]
    print(test)

    #residual feed
    result_file["residual_feed"] = pd.DataFrame()
    for i in range(len(nodes)):
        temp = [res_dom["feed"][i][t].X for t in time_steps]
        result_file["residual_feed"].insert(i, "Building_" + str(i+1), temp, True)

    #residual hp power electric
    result_file["residual_power_hp"] = pd.DataFrame()
    for i in range(len(nodes)):
        temp = [power_dom["HP"][i][t].X for t in time_steps]
        result_file["residual_power_hp"].insert(i, "Building_" + str(i), temp, True)

    #residual power pv
    result_file["residual_power_pv"] = pd.DataFrame()
    for i in range(len(nodes)):
        temp = [power_dom["PV"][i][t].X for t in time_steps]
        result_file["residual_power_pv"].insert(i, "Building_" + str(i), temp, True)

    #residual power bat
    result_file["residual_power_bat"] = pd.DataFrame()
    for i in range(len(nodes)):
        temp = [dch_dom["BAT"][i][t].X for t in time_steps]
        result_file["residual_power_bat"].insert(i, "Building_" + str(i), temp, True)

    #residual power bat
    result_file["residual_power_ev"] = pd.DataFrame()
    for i in range(len(nodes)):
        temp = [dch_dom["EV"][i][t].X for t in time_steps]
        result_file["residual_power_ev"].insert(i, "Building_" + str(i), temp, True)

    #residual soc TES
    result_file["residual_soc_TES"] = pd.DataFrame()
    for i in range(len(nodes)):
        temp = [soc_dom["TES"][i][t].X for t in time_steps]
        result_file["residual_soc_TES"].insert(i, "Building_" + str(i), temp, True)

    #residual soc BAT
    result_file["residual_soc_BAT"] = pd.DataFrame()
    for i in range(len(nodes)):
        temp = [soc_dom["BAT"][i][t].X for t in time_steps]
        result_file["residual_soc_BAT"].insert(i, "Building_" + str(i), temp, True)

    #trafo from grid
    result_file["from_grid"] = [power["from_grid"][t].X for t in time_steps]

    #trafo to grid
    result_file["to_grid"] = [power["to_grid"][t].X for t in time_steps]

    #sum residual feed for each building
    result_file ["residual_feed_sum"] = [sum(res_dom["feed"][i][t].X for t in time_steps)for i in range(len(nodes))]

    #sum residual power for each building
    result_file ["residual_power_sum"] = [sum(res_dom["power"][i][t].X for t in time_steps) for i in range(len(nodes))]

    #total feed in
    result_file["total_feed"] =sum(res_dom["feed"][i][t].X for t in time_steps for i in range(len(nodes)))

    #total power needed
    result_file["total_power"] =sum(res_dom["power"][i][t].X for t in time_steps for i in range(len(nodes)))

    store.save_dict(result_file, "results/Week"+str(idx)+"/central_optimization_week"+str(idx)+"_agg.npy")

    return result_file
