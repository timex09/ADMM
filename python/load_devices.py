#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 27 14:32:01 2020

@author: timowenzel, Sarah Henn
"""


def load_devices():

    devs = {}

    # BATTERY
    devs["bat"] = {}
    # Fronius 10.5 kWh, Symo Hybrid 4.0-3-S 9.6 # https://www.photovoltaik4all.de/fronius-energy-package-lithium-ionen-batterie-4,5-bis-12-kwh
    devs["bat"]["cat1"] = dict(cap=8.4, max_ch=5, max_dch=5, min_soc=0.1, max_soc=0.9, eta_bat=0.957)

    # Fronius 12.0 kWh, Symo Hybrid 4.0-3-S 9.6
    devs["bat"]["cat2"] = dict(cap=9.6, max_ch=6, max_dch=6, min_soc=0.1, max_soc=0.9, eta_bat=0.957)

    # Fronius 12.0 kWh, Symo Hybrid 4.0-3-S 9.6
    devs["bat"]["cat3"] = dict(cap=12, max_ch=9, max_dch=9, min_soc=0.1, max_soc=0.9, eta_bat=0.957)

    # PHOTOVOLTAICS
    # Honey Black TSM DD06M.05(II)
    # https://www.photovoltaik4all.de/media/pdf/3a/39/9d/Trina-HoneyBlack_DD06M_2020_A_Datasheet_HoneyBlackM_DD06M-05-II-_DE_01-2020.pdf
    # optical efficiency according to https: // www.homerenergy.com / products / pro / docs / 3.11 / solar_transmittance.html
    devs["PV"] = dict(eta_el_stc=0.199, area_real=1.6, eta_inv=0.96, t_cell_stc=25, G_stc=1000, t_cell_noct=44,
                      t_air_noct=20, G_noct=800, gamma=-0.003, eta_opt=0.9)

    # HEATPUMP
    devs["hp_air"] = {}
    # Viessmann Vitocal 200-S 6,1 kW
    # https://www.viessmann.de/de/wohngebaeude/waermepumpe/split-luft-wasser-waermepumpen/vitocal-200-s.html
    devs["hp_air"]["cat1"] = dict(cap=6.1, dT_max=15)

    # Viessmann Z016857 Luft-Wassser Wärmepumpe
    # https://www.meinhausshop.de/VIESSMANN-Z016857-Luft-Wasser-Waermepumpe-Vitocal-222-A-7-5-kW-Typ-AWOT-M-E-AC-221.A08-230-V
    devs["hp_air"]["cat2"] = dict(cap=7.5, dT_max=15)

    # Stiebel Luft/Wasser-Wärmepumpe WPL47ASR 26,46 kW, 232125
    # https://www.wolf-online-shop.de/Stiebel-Luft-Wasser-Waermepumpe-WPL47ASR-26-46-kW-232125::32855.html
    devs["hp_air"]["cat3"]  = dict(cap=26.46, dT_max=15)

    # ELECTRIC HEATER
    # Electric Heater 10000 W #2
    devs["eh"] = dict(cap=100)

    # THERMAL ENERGY STORAGE
    devs["tes"] = {}
    # Buderus Logalux SU300/5, 250 Liter
    devs["tes"]["cat1"] = dict(cap=0.25, dT_max=15, min_soc=0.1, eta_tes=0.98)

    # Buderus Logalux SU300/5, 300 Liter
    devs["tes"]["cat2"] = dict(cap=0.3, dT_max=15, min_soc=0.1, eta_tes=0.98)

    # Viessmann Vitocell 100-E, SVPA, 950 Liter
    devs["tes"]["cat3"] = dict(cap=0.95, dT_max=15, min_soc=0.1, eta_tes=0.995)

    # ELECTRIC VEHICLE
    # Renault Zoe
    devs["ev"] = dict(cap=35, eta_ch_ev=0.95, eta_dch_ev=0.97, min_soc=0.1, max_soc=0.9, max_ch_ev=45,
                      max_dch_ev=40)

    return devs

def map_devices(nodes, building_params, devs):

    for n in nodes:

        nodes[n]["PV"] = devs["PV"]
        nodes[n]["EV"] = devs["ev"]
        nodes[n]["EB"] = devs["eh"]


        if building_params["cat"][n] == 1:
            nodes[n]["BAT"] = devs["bat"]["cat1"]
            nodes[n]["HP"] = devs["hp_air"]["cat1"]
            nodes[n]["TES"] = devs["tes"]["cat1"]

        elif building_params["cat"][n] == 2:
            nodes[n]["BAT"] = devs["bat"]["cat2"]
            nodes[n]["HP"] = devs["hp_air"]["cat2"]
            nodes[n]["TES"] = devs["tes"]["cat2"]

        else:
            nodes[n]["BAT"] = devs["bat"]["cat3"]
            nodes[n]["HP"] = devs["hp_air"]["cat3"]
            nodes[n]["TES"] = devs["tes"]["cat3"]

    return nodes
