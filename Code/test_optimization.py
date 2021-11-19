"""
Den GridLineOptimizer testen
"""

from optimization import GridLineOptimizer as GLO
from battery_electric_vehicle import BatteryElectricVehicle as BEV
from household import Household as HH


resolution = 15
buses = 6
bevs = 5
bev_lst = list(range(bevs))
bus_lst = list(range(buses))
s_trafo = 150  #kVA

# BEVs
home_buses = [0, 1, 2, 3, 5]
start_socs = [20, 20, 30, 20, 40]
target_socs = [80, 70, 100, 90, 70]
target_times = [16, 16, 15, 18, 20]
bat_energies = [50, 50, 50, 50, 50]
bus_volts = [400-i/2 for i in bev_lst]

# Households
ann_dems = [3000, 35000, 3000, 4000, 3000, 3000]

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
    household.raise_demand(10, 17, 1000)
    household_list.append(household)

test = GLO(number_buses=buses, bevs=bev_list, resolution=resolution, s_trafo_kVA=s_trafo,
           households=household_list)

test.optimization_model.SOC.pprint()
test.optimization_model.household_currents.pprint()
test.display_max_current_constraint()
test.display_min_voltage_constraint()

# optimieren lassen
test.run_optimization_single_timestep(tee=True)

# Ergebnisse darstellen
test.plot_results()