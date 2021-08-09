#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 24 15:47:39 2020

@author: Sarah Henn
"""

# import python packages
import numpy as np
import pandas as pd
import python.clustering_medoid as cm
import python.pv_buildings as pv
import python.ev_param as ev_param
import demands.solar_power.solar as solar


def load_demands(options, devs):
    buildings = options["buildings"]
    datapoints = 8760
    path_input = options["path_input"]
    time_hor = options["time_hor"]
    tweeks = options["tweeks"]
    pv_share = options["pv_share"]
    ev_status = options["ev_status"]

    # Weather parameters
    weather_param = pd.read_csv(path_input + "weather.csv")
    weather_param["G_dif"] = 0.1 * weather_param["G_sol"]

    # Building parameters
    building_params = pd.read_csv(path_input + "Buildings_Bedburg.csv", delimiter=";")

    # Fill nodes dictionairy with demands
    demands = {}
    for i in range(buildings):
        demands[i] = {
            "elec_raw": np.loadtxt(
                open(path_input + "demands/elecLoadC" + str(building_params["names"][i]) + ".csv", "rb"),
                delimiter=",", skiprows=1, usecols=1) / 1000,  # kW, electricity demand
            "heat": np.loadtxt(
                open(path_input + "demands/demand_ClusterC" + str(building_params["names"][i]) + "_heatDemand.txt",
                     "rb"),
                delimiter=",", skiprows=3, usecols=1) / 1000,  # kW, heat demand
            "dhw": np.loadtxt(
                open(path_input + "demands/demand_ClusterC" + str(building_params["names"][i]) + "_dhwDemand.txt", "rb"),
                delimiter=",", skiprows=2, usecols=1) / 1000,  # kW, drinking hot water demand
        }

    # transform 15-min to hourly values
    elec = np.zeros(shape=(buildings, datapoints))
    for n in range(buildings):
        for t in range(datapoints):
            elec[n, t] = np.sum(demands[n]["elec_raw"][t * 4:(t * 4 + 4)]) / 4
        demands[n]["elec"] = elec[n, :]

    # get solar irradiance on tilted surface
    solar_irr = np.zeros(shape=(buildings, datapoints))
    for n in range(buildings):
        solar_irr[n] = solar.get_irrad_profile(weather_param, building_params["azimuth"][n],
                                               building_params["elevation"][n])
        demands[n]["solar_irr"] = solar_irr[n]

    # calculate time variant efficiency
    t_cell = np.zeros(shape=(datapoints, 1))
    eta_el = np.zeros(shape=(datapoints, 1))
    for t in range(datapoints):
        t_cell[t] = (weather_param["T_air"][t] + (devs["PV"]["t_cell_noct"] - devs["PV"]["t_air_noct"])
                     * (weather_param["G_sol"][t] / devs["PV"]["G_noct"])
                     * (1 - (devs["PV"]["eta_el_stc"] * (1 - devs["PV"]["gamma"] * devs["PV"]["t_cell_stc"])) /
                        devs["PV"]["eta_opt"])) / \
                    ((1 + (devs["PV"]["t_cell_noct"] - devs["PV"]["t_air_noct"])
                      * (weather_param["G_sol"][t] / devs["PV"]["G_noct"]) *
                      ((devs["PV"]["gamma"] * devs["PV"]["eta_el_stc"]) / devs["PV"]["eta_opt"])))

        eta_el[t] = devs["PV"]["eta_el_stc"] * (1 + devs["PV"]["gamma"] * (t_cell[t] - devs["PV"]["t_cell_stc"]))

        devs["PV"]["eta_el"] = eta_el[:]

    # calculate PV power
    pv_power = np.zeros(shape=(buildings, datapoints))
    for n in range(buildings):
        for t in range(datapoints):
            pv_power[n, t] = demands[n]["solar_irr"][t] * devs["PV"]["eta_el"][t] * devs["PV"]["area_real"] \
                             * building_params["modules"][n] / 1000  # kW,       PV output

        demands[n]["pv_power"] = pv_power[n, :]

    # calculate availability of pv system considerung pv_share
    pv_avail = pv.pv_share(pv_share, buildings)

    # adjust time steps for typical weeks: 52*7*24 = 8736 instead of 8760
    demands_adj = {
        "T_air": weather_param["T_air"][0:-24]}
    for n in range(buildings):
        demands_adj[n] = {"elec": demands[n]["elec"][0:-24],
                          "heat": demands[n]["heat"][0:-24],
                          "dhw": demands[n]["dhw"][0:-24],
                          "pv_power": demands[n]["pv_power"][0:-24]}

    # Prepare clustering
    inputs_clustering = []

    for n in range(buildings):
        inputs_clustering.append(demands_adj[n]["elec"])
        inputs_clustering.append(demands_adj[n]["heat"])
        inputs_clustering.append(demands_adj[n]["dhw"])
        inputs_clustering.append(demands_adj[n]["pv_power"])

    inputs_clustering.append(demands_adj["T_air"])

    # Perform clustering
    (inputs, nc, z, inputs_transformed) = cm.cluster(np.array(inputs_clustering),
                                                     number_clusters=tweeks, len_day=time_hor)

    # Fill nodes dictionairy with clustered demand an weather data and pv panels data
    # Todo: Noch nicht die optimale Lösung PV mitzuclustern und dann gleich 0 zu setzen... Wie geht das besser?
    nodes = {}
    for n in range(buildings):
        nodes[n] = {
            "elec": inputs[4 * n].T,
            "heat": inputs[4 * n + 1].T,
            "dhw": inputs[4 * n + 2].T,
            "pv_power": pv_avail[n] * inputs[4 * n + 3].T,
            "T_air": inputs[-1].T,
            "weights": nc
        }

    # Check for small demand values
    for n in range(buildings):
        for t, m in zip(range(time_hor), range(tweeks)):
            if nodes[n]["elec"][t, m] < 0.01:
                nodes[n]["elec"][t, m] = 0
            if nodes[n]["heat"][t, m] < 0.01:
                nodes[n]["heat"][t, m] = 0
            if nodes[n]["dhw"][t, m] < 0.01:
                nodes[n]["dhw"][t, m] = 0

    # Calculation of Coefficient of Power with temperatures of first building
    for n in nodes:
        nodes[n]["COP_sh"] = np.zeros_like(nodes[n]["heat"])
        nodes[n]["COP_dhw"] = np.zeros_like(nodes[n]["dhw"])
        for m in range(tweeks):
            for t in range(time_hor):
                nodes[n]["COP_sh"][t, m] = 0.4 * (273.15 + building_params["T_sup_sh"][n]) / \
                                           (building_params["T_sup_sh"][n] - nodes[n]["T_air"][t, m])
                nodes[n]["COP_dhw"][t, m] = 0.4 * (273.15 + building_params["T_sup_dhw"][n]) / \
                                            (building_params["T_sup_dhw"][n] - nodes[n]["T_air"][t, m])

    # Todo: ev_param Script checken. Auf den ersten Blick haben dem_arrive und dem_leave nicht die gleiche Datenbasis
    # Electric Vehicles
    counter = 0
    ev_data = ev_param.ev_load(path_input)

    ev_avail = np.zeros(shape=(time_hor, 1))
    ev_dem_arrive = np.zeros(shape=(time_hor, 1))
    ev_dem_leave = np.zeros(shape=(time_hor, 1))
    for n in nodes:
        if ev_status and building_params["ev_exists"][n] == 1:
            for t in range(time_hor):
                ev_avail[t] = ev_data["avail"][counter, t]
                ev_dem_arrive[t] = ev_data["dem_arrive"][t, counter]
                ev_dem_leave[t] = ev_data["dem_leave"][t, counter]
            nodes[n]["ev_avail"] = ev_avail[:]
            nodes[n]["ev_dem_arrive"] = ev_dem_arrive[:]
            nodes[n]["ev_dem_leave"] = ev_dem_leave[:]
            counter = counter + 1

        else:
            nodes[n]["ev_avail"] = ev_avail
            nodes[n]["ev_dem_arrive"] = ev_dem_arrive
            nodes[n]["ev_dem_leave"] = ev_dem_leave

    return nodes, building_params
