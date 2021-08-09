#!/usr/bin/env python
# coding=utf-8
"""
Example script on how to generate a stochastic electric load profile
"""

# #  Seed is for testing purpose, only!
import random as rd
# rd.seed(1)

import copy
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

import richardsonpy.classes.occupancy as occ
import richardsonpy.functions.change_resolution as cr
import richardsonpy.functions.load_radiation as loadrad
import richardsonpy.classes.electric_load as eload


def example_stoch_el_load_sfh(nb_occ,nb, do_plot=False):
    #  Total number of occupants in apartment
    nb_occ = nb_occ

    timestep = 900  # in seconds

    #  Generate occupancy object (necessary as input for electric load gen.)
    occ_obj = occ.Occupancy(number_occupants=nb_occ)

    #  Get radiation (necessary for lighting usage calculation)
    (q_direct, q_diffuse) = loadrad.get_rad_from_try_path()

    #  Convert 3600 s timestep to given timestep
    q_direct = cr.change_resolution(q_direct, old_res=3600, new_res=timestep)
    q_diffuse = cr.change_resolution(q_diffuse, old_res=3600, new_res=timestep)

    #  Generate stochastic electric power object
    el_load_obj = eload.ElectricLoad(occ_profile=occ_obj.occupancy,
                                     total_nb_occ=nb_occ,
                                     q_direct=q_direct,
                                     q_diffuse=q_diffuse,
                                     timestep=timestep)

    #  Copy occupancy object, before changing its resolution
    #  (occ_obj.occupancy is the pointer to the occupancy profile array)
    occ_profile_copy = copy.copy(occ_obj.occupancy)

    #  Change resolution of occupancy object (to match
    #  resolution of el. load profile; necessary for plotting)
    occ_profile_copy = cr.change_resolution(values=occ_profile_copy,
                                            old_res=600,
                                            new_res=timestep)

    #  Calculate el. energy in kWh
    energy_el_kwh = sum(el_load_obj.loadcurve) * timestep / (3600 * 1000)

    print('Electric energy demand in kWh: ')
    print(energy_el_kwh)
    
    print(el_load_obj.loadcurve)
    
    loadProfile = pd.DataFrame(el_load_obj.loadcurve)
    loadProfile.index = loadProfile.index * timestep
    
    #loadProfile.to_csv("sfh.csv", header = 'SFH' + str(nb_occ))
    loadProfile.to_csv("sfh"+str(nb)+"_" +str(nb_occ) +"persons.csv", index_label='time in s', header = ['power in W'])

    # if do_plot:
    #     #  Generate time array for plotting
    #     timesteps = int(8760 * 3600 / timestep)  # Number of timesteps per day

    #     time_array = np.arange(0, timesteps * timestep, timestep) / 3600

    #     fig = plt.figure()
    #     fig.add_subplot(211)
    #     plt.plot(time_array, occ_profile_copy[0:timesteps])
    #     plt.xlabel('Timestep in hours')
    #     plt.ylabel('Number of active occupants')
    #     plt.xlim([0, 8760])

    #     fig.add_subplot(212)
    #     plt.plot(time_array, el_load_obj.loadcurve[0:timesteps])
    #     plt.xlabel('Timestep in hours')
    #     plt.ylabel('Electric power in W')
    #     plt.xlim([0, 8760])

    #     plt.tight_layout()
    #     plt.show()
    #     #plt.savefig('DiversityIndex.pdf', dpi=1200)
    #     #plt.savefig('DiversityIndex.png', dpi=1200, transparent=True)
    #     plt.close()


def example_stoch_el_load_mfh_flat(nb_occ):
    #  Total number of occupants in apartment
    nb_occ = nb_occ

    timestep = 900  # in seconds

    #  Generate occupancy object (necessary as input for electric load gen.)
    occ_obj = occ.Occupancy(number_occupants=nb_occ)

    #  Get radiation (necessary for lighting usage calculation)
    (q_direct, q_diffuse) = loadrad.get_rad_from_try_path()

    #  Convert 3600 s timestep to given timestep
    q_direct = cr.change_resolution(q_direct, old_res=3600, new_res=timestep)
    q_diffuse = cr.change_resolution(q_diffuse, old_res=3600, new_res=timestep)

    #  Generate stochastic electric power object
    el_load_obj = eload.ElectricLoad(occ_profile=occ_obj.occupancy,
                                     total_nb_occ=nb_occ,
                                     q_direct=q_direct,
                                     q_diffuse=q_diffuse,
                                     timestep=timestep,
                                     is_sfh=False)

    #  Copy occupancy object, before changing its resolution
    #  (occ_obj.occupancy is the pointer to the occupancy profile array)
    occ_profile_copy = copy.copy(occ_obj.occupancy)

    #  Change resolution of occupancy object (to match
    #  resolution of el. load profile; necessary for plotting)
    occ_profile_copy = cr.change_resolution(values=occ_profile_copy,
                                            old_res=600,
                                            new_res=timestep)

    #  Calculate el. energy in kWh
    energy_el_kwh = sum(el_load_obj.loadcurve) * timestep / (3600 * 1000)

    print('Electric energy demand in kWh: ')
    print(energy_el_kwh)
    
    print(el_load_obj.loadcurve)
    
    loadProfile = pd.DataFrame(el_load_obj.loadcurve)
    loadProfile.index = loadProfile.index * timestep
    
    return(loadProfile)

if __name__ == '__main__':
    nb_sfh = 5
    nb_mfh = 2
    nb_flats = 4
    
    
    for i in range(1,nb_sfh+1):
        nb_occ =  rd.randrange(2, 5, 1)
        print(nb_occ)
        example_stoch_el_load_sfh(nb_occ,nb=i, do_plot=False)
    
    sum_mfh = pd.DataFrame() 
    nb_occ_total = 0

    for i in range(1, nb_mfh+1):
        for j in range (1, nb_flats+1):    
            nb_occ =  rd.randrange(2, 5, 1)
            el_load_flat = example_stoch_el_load_mfh_flat(nb_occ)
            sum_mfh = sum_mfh.add(el_load_flat, fill_value=0)
            nb_occ_total += nb_occ
        sum_mfh.to_csv("mfh"+str(i)+"_" +str(nb_occ_total) +"persons.csv", index_label='time in s', header = ['power in W'])
            
            
            
            
            
            
            
            
            
            