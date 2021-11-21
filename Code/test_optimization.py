"""
Den GridLineOptimizer testen
"""

from optimization import GridLineOptimizer as GLO
from battery_electric_vehicle import BatteryElectricVehicle as BEV
from household import Household as HH


resolution = 60
buses = 10
bevs = 10
bev_lst = list(range(bevs))
bus_lst = list(range(buses))
s_trafo = 150  #kVA

# BEVs
home_buses = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
start_socs = [20, 20, 30, 20, 40, 25, 30, 20, 30, 20]
target_socs = [80, 70, 100, 90, 70, 85, 90, 80, 90, 80]
target_times = [16, 16, 15, 18, 20, 12, 18, 19, 22, 20]
bat_energies = [50, 50, 50, 50, 50, 50, 50, 50, 50, 50]
bus_volts = [400-i/2 for i in bev_lst]

# Households
ann_dems = [3000, 3500, 3000, 4000, 3000, 3000, 3500, 4000, 3000, 3000]

# BEVs erzeugen
bev_list = []
for car in bev_lst:
    bev = BEV(soc_start=start_socs[car], soc_target=target_socs[car],
              t_target=target_times[car], e_bat=bat_energies[car],
              resolution=resolution, bus_voltage=bus_volts[car],
              home_bus=home_buses[car])
    bev_list.append(bev)

# Households erzeugen
household_list = []
for bus in bus_lst:
    household = HH(home_bus=bus, annual_demand=ann_dems[bus], resolution=resolution)
    #household.raise_demand(12, 21, 7500)
    household_list.append(household)

test = GLO(number_buses=buses, bevs=bev_list, resolution=resolution, s_trafo_kVA=s_trafo,
           households=household_list)


# optimieren lassen
test.run_optimization_single_timestep(tee=True)

test.optimization_model.SOC.pprint()
#test._prepare_next_timestep()

# Ergebnisse darstellen
test.plot_results()
