"""
Author: André Ulrich
--------------------
Den GridLineOptimizer testen

Versionsgeschichte:
V.1: jetzt kann unterschieden werden, ob man mit oder ohne rolling horizon optimieren möchte

V.2:
"""
import random

from optimization import GridLineOptimizer as GLO
from battery_electric_vehicle import BatteryElectricVehicle as BEV
from household import Household as HH
from EMO import *

import matplotlib.pyplot as plt

ROLLING = True#'experimental'
random_wishes = False
use_emo = False # run EMO simulation to verify the optimization results
emo_unoptimized = False # run EMO sinulation without optimization (BEVs charge according to P(SOC) curve) but P(U) controling
emo_uncontrolled = False # run EMO simulation without optimization and without controlling

#========================================================
# define scenario
#========================================================

seed = 5 # for creating reproducible "random" numbers
resolution = 15 # resolution in minutes
horizon = 24 # time horizon [h]
buses = 6 # buses on the grid line (excluding trafo lv and mv and slack)
bevs = 2 # buses with charging station (makes no sense to choose greater than buses)
p_trafo = 15  # power of transformer [kVA]
bev_lst = list(range(bevs)) # for iterating purposes
bus_lst = list(range(buses)) # for iterating purposes

#==== if random generated wishes, these are the characteristics of the random distribution
# you need to take care that these make sense (e.g. target soc is greater than start soc)
# as software won't raise exception =========================================================
# soc [%]
start_socs_mean = 10
start_socs_deviation_plus = 0
start_socs_deviation_minus = 0

target_socs_mean = 100
target_socs_deviation_plus = 0
target_socs_deviation_minus = 0

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
    home_buses = [i for i in range(bevs)]
    start_socs = [start_socs_mean + random.randint(-1 * start_socs_deviation_minus, start_socs_deviation_plus) for _ in range(bevs)]
    target_socs = [target_socs_mean + random.randint(-1 * target_socs_deviation_minus, target_socs_deviation_plus) for _ in range(bevs)]
    target_times = [target_times_mean + random.randint(-1 * target_times_deviation_minus, target_times_deviation_plus) for _ in range(bevs)]
    start_times = [start_times_mean + random.randint(-1 * start_times_deviation_minus, start_times_deviation_plus) for _ in range(bevs)]
    bat_energies = [50 for _ in range(bevs)]

else: # create them on your own (length of list must equal bevs)
    home_buses = [0, 5]
    start_socs = [50, 20]
    target_socs = [100, 100]
    target_times = [20, 20]
    start_times = [16, 16]
    bat_energies = [50, 50]


#==== Households annual demand [kWh] ====================================================================
ann_dems = [3500 for _ in range(buses)]

#===========================================================================
# create objects for optimization
#===========================================================================

#==== create the BEV instances ==============================================================
bev_list = []
for car in bev_lst:
    bev = BEV(soc_start=start_socs[car], soc_target=target_socs[car],
              t_target=target_times[car], e_bat=bat_energies[car],
              resolution=resolution, home_bus=home_buses[car],
              t_start=start_times[car])
    bev_list.append(bev)

#==== create the Households instances =======================================================
household_list = []
for bus in bus_lst:
    household = HH(home_bus=bus, annual_demand=ann_dems[bus], resolution=resolution)
    #household.raise_demand(11, 19, 23800) # raise demand to simulate additional electric loads if you like
    household_list.append(household)

#==== choose additional setup for optimizer =================================================
#GLO.set_options('log results', True)
#GLO.set_options('fairness', 2)
#GLO.set_options('equal SOCs', 0.1)
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
    test.run_optimization_single_timestep(tee=True)
    test.optimization_model.SOC.pprint()
    test.plot_all_results(marker=None, save=False, usetex=True, compact_x=True, export_data=True)
    #test.plot_I_results(marker=None, save=True, usetex=True, compact_x=True)
    #test.plot_SOC_results(marker=None, save=True, usetex=True, compact_x=True)
    test.export_socs_fullfillment()

else:
    test.run_optimization_rolling_horizon(tee=False, complete_horizon=24)
    test.plot_all_results(marker=None, save=False, usetex=True, compact_x=True, export_data=True)
    test.export_socs_fullfillment(rolling=True)

if use_emo: # falls Rolling, dann sind ja schon die SOCs der BEVs angepasst woreden => müssen wieder auf SOC_start gesetzt werden
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
        sim_handler_1.run_GLO_sim(hh_data, wb_data, int(horizon * 60 / resolution), parallel=False)

    else:
        if not emo_uncontrolled:
            sim_handler_1.run_unoptimized_sim(hh_data, bev_list, int(horizon * 60 / resolution), control=True)
            test.export_socs_fullfillment(optimized=False)

        else:
            sim_handler_1.run_unoptimized_sim(hh_data, bev_list, int(24 * 60 / resolution), control=False)
            test.export_socs_fullfillment(optimized=False)

#==== Visualize the simulation results ==============================================
    sim_handler_1.plot_EMO_sim_results(freq=resolution, element='buses', legend=False, marker=None,
                                       save=True, usetex=True, compact_x=True)
    sim_handler_1.plot_EMO_sim_results(freq=resolution, element='lines', legend=False, marker=None,
                                       save=True, usetex=True, compact_x=True)
    sim_handler_1.plot_EMO_sim_results(freq=resolution, element='trafo', legend=False, marker=None,
                                       save=True, usetex=True, compact_x=True)

    sim_handler_1.export_sim_results('trafo', res_min=resolution)
    sim_handler_1.export_sim_results('buses', res_min=resolution)


# if ROLLING == True:
#     test.run_optimization_rolling_horizon(24, tee=False)
#     for key in test.results_I:
#         print(test.results_I[key])
#
#     for i in home_buses:#range(len(bev_lst)):
#         plt.plot(range(len(test.results_I[i])), test.results_SOC[i])#, marker='o')
#     plt.show()
#
#     for i in home_buses:#range(len(bev_lst)):
#         plt.plot(range(len(test.results_I[i])), test.results_I[i], label=f'Current to BEV at node {i}')#, marker='o')
#     plt.legend()
#     plt.show()
#
#     for bev in bev_list:
#         print(f'Verlauf der SOCs des BEV an Knoten {bev.home_bus}', bev.soc_list, '\n')
#
#
# if ROLLING == 'experimental':
#     res_ges = []
#     for t in range(int(24*60/resolution)):
#         print(t)
#         test.run_optimization_single_timestep(tee=False)
#         test._store_results()
#         res = test.export_I_results()
#         res_ges.append(res)
#         test._prepare_next_timestep()
#         test._setup_model()
#
#     for item in res_ges:
#         for key in item:
#             print(item[key])
#         print('############################################')



