#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 29 15:03:17 2020

@author: Leo
"""

import numpy as np
import pandas as pd
import math


def save_dict(dictionary, file_code):
    
    np.save(file_code, dictionary)

    return



def read_dict(file_code):
    
    read_dictionary = np.load(file_code,allow_pickle='TRUE').item()
    
    return read_dictionary


def dyn_net_factor(powerLine):
        
    factor = powerLine*1.5
    
    return factor



"""
dictionary_WirtGes = read_dict("result_dict.npy")

dictionary_WirtGes_zen = read_dict("result_dict_zenBat.npy")

dictionary_NetSingle = read_dict("result_dict_net_single.npy")

dictionary_NetSumNoFac = read_dict("result_dict_net_all_without_factor.npy")

dictionary_NetSumLinFac = read_dict("result_dict_net_all.npy")


df = pd.DataFrame()

df["WirtGes"] = dictionary_WirtGes["PowerLines"][(1,212)]

df["NetSingle"] = dictionary_NetSingle["PowerLines"][(1,212)]

df["NetSumNoFac"] = dictionary_NetSumNoFac["PowerLines"][(1,212)]

df["NetSumLinFac"] = dictionary_NetSumLinFac["PowerLines"][(1,212)]


boxplot = df.boxplot()


df_temp = pd.DataFrame()

df_temp["Residual_power"] = dictionary_WirtGes["residual_df_power"]["residual_power"]
df_temp["Residual_feed"] = dictionary_WirtGes["residual_df_feed"]["residual_feed"]
df_temp["from_grid"] = np.asarray(dictionary_WirtGes["power_from_grid"])*1000


df_temp["Residual_power_zen"] = dictionary_WirtGes_zen["residual_df_power"]["residual_power"]
df_temp["Residual_feed_zen"] = dictionary_WirtGes_zen["residual_df_feed"]["residual_feed"]
df_temp["from_grid_zen"] = np.asarray(dictionary_WirtGes_zen["power_from_grid"])

df_temp.to_excel("Energieaustausch_Vergleich.xlsx")

df = pd.DataFrame()
df["dezentral"] = dictionary_WirtGes["PowerLines"][(1,212)]
df["zentral"] = dictionary_WirtGes_zen["PowerLines"][(1,212)]

boxplot = df.boxplot()
"""