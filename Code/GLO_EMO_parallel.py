"""
Author: André Ulrich
--------------------
Test parallelization of optimization (GridLineOptimizer) and grid simulation (EMO).
Optimization is running in own separate process and sends the results of optimized
BEV charging currents (only the ones of the current timestep) through a queue
to the EMO where an according grid silation is run.

Version history (only the most relevant points, full history is available on github):
-------------------------------------------------------------------------------------------------
V.1: first working parallelization approach

V.2: separate thread for reading new results from queue

all the other commits in much more detail are available here:
https://github.com/AndreGismo/Masterarbeit/tree/submission)
"""
import os

from EMO import *
from optimization import GridLineOptimizer as GLO
from battery_electric_vehicle import BatteryElectricVehicle as BEV
from household import Household as HH

from matplotlib import pyplot as plt
import matplotlib.pyplot as plt
import multiprocessing as mp
import threading as thr
import os
import time

CPUS = os.cpu_count()

#### GridLineOptimizer ########################################################
resolution = 15
buses = 6
bevs = 6
bev_lst = list(range(bevs))
bus_lst = list(range(buses))
s_trafo = 150  #kVA

t_steps = 96

t_counter = 0

# BEVs
home_buses = [0, 1, 2, 3, 4, 5]
start_socs = [20, 20, 30, 20, 25, 40]
target_socs = [80, 70, 80, 90, 80, 70]
target_times = [10, 18, 18, 18, 18, 18]
start_times = [2, 2, 2, 2, 2, 2]
bat_energies = [50, 50, 50, 50, 50, 50]

# Households
ann_dems = [3000, 3500, 3000, 4000, 3000, 3000]

# BEVs erzeugen
bev_list = []
for car in bev_lst:
    bev = BEV(soc_start=start_socs[car], soc_target=target_socs[car],
              t_target=target_times[car], e_bat=bat_energies[car],
              resolution=resolution, home_bus=home_buses[car],
              t_start=start_times[car])
    bev_list.append(bev)

# Households erzeugen
household_list = []
for bus in bus_lst:
    household = HH(home_bus=bus, annual_demand=ann_dems[bus], resolution=resolution)
    #household.raise_demand(11, 19, 23500)
    household_list.append(household)

test = GLO(number_buses=buses, bevs=bev_list, resolution=resolution, trafo_power=s_trafo,
           households=household_list, horizon_width=24)

# export grid as excel
grid_excel_file = 'optimized_grid'
test.export_grid(grid_excel_file)
grid_specs = test.get_grid_specs()
hh_data = test.export_household_profiles()
#wb_data = test.export_I_results() erst nach Optimierung sinnvoll! sonst nur None


#### EMO ######################################################################
system_1 = Low_Voltage_System(line_type='NAYY 4x120 SE', transformer_type="0.25 MVA 10/0.4 kV")
system_1.grid_from_GLO('optimized_grid.xlsx', grid_specs)

sim_handler_1 = Simulation_Handler(system_1,
                                    start_minute=60 * 12,
                                    end_minute=60 * 12 + 24 * 60,
                                    rapid=False)

queue = mp.Queue()

#func = test.run_optimization_single_timestep
def func_opt(tee, marker, queue):
    pid = os.getpid()
    global t_counter
    for t in range(t_steps):
        print('optimization round', t, 'running in process', pid)
        test.run_optimization_fixed_horizon(tee=tee)
        I_res = test.export_current_I_results(1)
        queue.put(I_res)
        test._store_results()
        test._prepare_next_timestep()
        test._setup_model()
        time.sleep(0.75)

    queue.put('done')


def func_sim(queue):
    pid = os.getpid()
    global t_counter

    # function to be run in own thread, to get the data (I_res) out of the
    # queue
    def monitor_queue():
        nonlocal queue, res_I
        while True:
            res_I = queue.get()

    # first run optimization for first timestep to have results
    test.run_optimization_fixed_horizon(tee=False)
    res_I = test.export_current_I_results(1)

    last_I_res = res_I
    # start the thread to monitor the queue
    thr.Thread(target=monitor_queue, daemon=True).start()

    while True:
        if res_I == last_I_res:
            print('no new results, will continue with the old one')

        else:
            print('received new results!')
            last_I_res = res_I

        if res_I == 'done':
            break

        # run simulation with only the results for the first timestep
        sim_handler_1.run_GLO_sim(hh_data, res_I, parallel=True)
        time.sleep(0.75)


if __name__ == '__main__':
    # create and start the process
    p_opt = mp.Process(target=func_opt, kwargs={'tee': False, 'marker': 'x', 'queue': queue}, daemon=True)
    p_opt.start()
    # .. and start parallel execution
    func_sim(queue)

    p_opt.join()
    print('done!')

    for i in range(6):
        plt.plot(range(len(sim_handler_1.res_GLO_sim_trafo)), sim_handler_1.res_GLO_sim_U[i])
    plt.show()

    plt.plot(range(len(sim_handler_1.res_GLO_sim_trafo)), sim_handler_1.res_GLO_sim_trafo)
    plt.show()

    test.plot_all_results(marker=None, save=False, usetex=False, compact_x=True, export_data=True)
