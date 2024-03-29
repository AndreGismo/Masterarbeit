"""
Author: André Ulrich
--------------------
testcase for optimization with GridLineOptimizer and validation of results
with grid simulation with EMO.

Uses all the functionalities of GridLineOptimizer and EMOs Simulation_Handler.

Version history (only the most relevant points, full history is available on github):
-------------------------------------------------------------------------------------------------
V.1: first working test case formulated

V.2: some switches at beginning to select what exactly to do

all the other commits in much more detail are available here:
https://github.com/AndreGismo/Masterarbeit/tree/submission)
"""
import random

from optimization import GridLineOptimizer as GLO
from battery_electric_vehicle import BatteryElectricVehicle as BEV
from household import Household as HH
from EMO import *

import matplotlib.pyplot as plt
import numpy as np

ROLLING = False # use rolling horizon instead of one fixed horizon for optimization
random_wishes = True # random generated customer wishes (usefull, when using lost of BEVs)
use_emo = True # run EMO simulation to verify the optimization results
emo_unoptimized = False # run EMO sinulation without optimization (BEVs charge according to P(SOC) curve) but P(U) controling
emo_uncontrolled = False # run EMO simulation without optimization and without controlling


#========================================================
# define scenario
#========================================================

seed = 7 # for creating reproducible "random" numbers
resolution = 15 # resolution in minutes
horizon = 24 # time horizon [h]
buses = 40 # buses on the grid line (excluding trafo lv and mv and slack)
bevs = 40 # buses with charging station (makes no sense to choose greater than buses)
p_trafo = 250  # power of transformer [kVA]
bev_lst = list(range(bevs)) # for iterating purposes
bus_lst = list(range(buses)) # for iterating purposes

#==== if random generated wishes, these are the characteristics of the random distribution
# you need to take care that these make sense (e.g. target soc is greater than start soc)
# as software won't raise exception =========================================================
# soc [%]
start_socs_mean = 30
start_socs_deviation_plus = 10
start_socs_deviation_minus = 10

target_socs_mean = 80
target_socs_deviation_plus = 20
target_socs_deviation_minus = 20

# time [h]
start_times_mean = 10
start_times_deviation_plus = 2
start_times_deviation_minus = 2

target_times_mean = 19
target_times_deviation_plus = 4
target_times_deviation_minus = 4


#==== BEVs customer wishes ==============================================================================
if random_wishes:
    random.seed(seed)
    np.random.seed(seed)
    home_buses = np.random.choice([i for i in range(buses)], size=bevs, replace=False)
    start_socs = [start_socs_mean + random.randint(-1 * start_socs_deviation_minus, start_socs_deviation_plus) for _ in range(bevs)]
    target_socs = [target_socs_mean + random.randint(-1 * target_socs_deviation_minus, target_socs_deviation_plus) for _ in range(bevs)]
    target_times = [target_times_mean + random.randint(-1 * target_times_deviation_minus, target_times_deviation_plus) for _ in range(bevs)]
    start_times = [start_times_mean + random.randint(-1 * start_times_deviation_minus, start_times_deviation_plus) for _ in range(bevs)]
    bat_energies = [50 for _ in range(bevs)]

else: # create them on your own (length of list must equal bevs)
    home_buses = [0, 5]
    start_socs = [10, 30]
    target_socs = [80, 80]
    target_times = [21, 21]
    start_times = [12, 16]
    bat_energies = [50, 50]


#==== Households annual demand [kWh] ====================================================================
ann_dems = [3500 for _ in range(buses)]

#===========================================================================
# create objects for optimization
#===========================================================================

#==== create the BEV instances ==============================================================
bev_list = []
for car in range(len(home_buses)):#home_buses:
    bev = BEV(soc_start=start_socs[car], soc_target=target_socs[car],
              t_target=target_times[car], e_bat=bat_energies[car],
              resolution=resolution, home_bus=home_buses[car],
              t_start=start_times[car])
    bev_list.append(bev)

#==== create the Households instances =======================================================
household_list = []
for bus in bus_lst:
    household = HH(home_bus=bus, annual_demand=ann_dems[bus], resolution=resolution)
    household.raise_demand(8, 16, 1700) # raise demand to simulate additional electric loads if you like
    household_list.append(household)

#==== choose additional setup for optimizer =================================================
#GLO.set_options('log results', True)
#GLO.set_options('fairness', 0)
GLO.set_options('equal SOCs', 0)
#GLO.set_options('equal products', True)
#GLO.set_options('atillas constraint', True)

#==== create optimizer instance =============================================================
test = GLO(number_buses=buses, bevs=bev_list, resolution=resolution, trafo_power=p_trafo,
           households=household_list, horizon_width=horizon)

#============================================================================================
# run the actual optimization
#============================================================================================

#==== standalone optimization ==============================================================
if not ROLLING:
    test.run_optimization_fixed_horizon(tee=True)
    test.optimization_model.SOC.pprint()
    # print results, use export_date=True to export results
    test.plot_all_results(marker=None, save=False, usetex=False, compact_x=True, export_data=True)
    #test.plot_I_results(marker=None, save=True, usetex=True, compact_x=True)
    #test.plot_SOC_results(marker=None, save=True, usetex=True, compact_x=True)
    test.export_socs_fullfillment()

else:
    test.run_optimization_rolling_horizon(tee=False, complete_horizon=horizon)
    test.plot_all_results(marker=None, save=False, usetex=False, compact_x=True, export_data=True)
    test.export_socs_fullfillment()

if use_emo:
#==== optimization + validation of results by using emo simulation: first prepare the data
# of the optimizer to be communicated to the emo-objects =====================================
    grid_excel_file = 'optimized_grid'
    test.export_grid(grid_excel_file)
    grid_specs = test.get_grid_specs()
    hh_data = test.export_household_profiles()
    wb_data = test.export_I_results()

#==== create the emo grid object and pass it the data of the grid that was used by the optimizer
    system_1 = Low_Voltage_System(line_type='NAYY 4x120 SE', transformer_type="0.25 MVA 10/0.4 kV")
    system_1.grid_from_GLO('optimized_grid.xlsx', grid_specs)

    sim_handler_1 = Simulation_Handler(system_1, start_minute=60 * 12, end_minute=60 * 12 + 24 * 60, rapid=False)

#==== start the emo net simulation =============================================================
    if not emo_unoptimized:
        sim_handler_1.run_GLO_sim(hh_data, wb_data, parallel=False)#, int(horizon * 60 / resolution), parallel=False)
        print('Starte Netzsimulation mit Optimierungsergebnissen')

    else:
        if not emo_uncontrolled:
            sim_handler_1.run_unoptimized_sim(hh_data, bev_list, int(horizon * 60 / resolution), control=True)
            test.export_socs_fullfillment(optimized=False)
            print('Starte Netzsimulation ohne Optimierungsergebnissen, nutze P(U)-Regelung')

        else:
            sim_handler_1.run_unoptimized_sim(hh_data, bev_list, int(horizon * 60 / resolution), control=False)
            test.export_socs_fullfillment(optimized=False)
            print('Starte Netzsimulation ohne Optimierungsergebnissen, freies Laden')

#==== Visualize the simulation results ==============================================
    sim_handler_1.plot_EMO_sim_results(freq=resolution, element='buses', legend=False, marker=None,
                                       save=True, usetex=True, compact_x=True)
    sim_handler_1.plot_EMO_sim_results(freq=resolution, element='lines', legend=False, marker=None,
                                       save=True, usetex=True, compact_x=True)
    sim_handler_1.plot_EMO_sim_results(freq=resolution, element='trafo', legend=False, marker=None,
                                       save=True, usetex=True, compact_x=True)

    sim_handler_1.export_sim_results('trafo', res_min=resolution)
    sim_handler_1.export_sim_results('buses', res_min=resolution)



