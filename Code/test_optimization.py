"""
Author: André Ulrich
--------------------
Den GridLineOptimizer testen

Versionsgeschichte:
V.1: jetzt kann unterschieden werden, ob man mit oder ohne rolling horizon optimieren möchte

V.2:
"""

from optimization import GridLineOptimizer as GLO
from battery_electric_vehicle import BatteryElectricVehicle as BEV
from household import Household as HH

import matplotlib.pyplot as plt

ROLLING = False

resolution = 15
buses = 6
bevs = 6
bev_lst = list(range(bevs))
bus_lst = list(range(buses))
s_trafo = 50  #kVA

# BEVs
home_buses = [0, 1, 2, 3, 4, 5]
start_socs = [20, 20, 30, 20, 25, 40]
target_socs = [80, 70, 80, 90, 80, 70]
target_times = [10, 16, 18, 18, 17, 20]
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

test = GLO(number_buses=buses, bevs=bev_list, resolution=resolution, s_trafo_kVA=s_trafo,
           households=household_list, horizon_width=24)


# optimieren lassen
if not ROLLING:
    test.run_optimization_single_timestep(tee=True)
    test.optimization_model.SOC.pprint()
    test.plot_results(marker='o')
    test.optimization_model.occupancy_times.pprint()

else:
    test.run_optimization_rolling_horizon(24, tee=False)
    for key in test.results_I:
        print(test.results_I[key])

    for i in range(len(bev_lst)):
        plt.plot(range(len(test.results_I[0])), test.results_SOC[i], marker='o')
    plt.show()



