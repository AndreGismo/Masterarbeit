"""
Den GridLineOptimizer testen
"""

from optimization import GridLineOptimizer as GLO
from battery_electric_vehicle import BatteryElectricVehicle as BEV

resolution = 30
buses = 6
bus_lst = list(range(buses))
s_trafo = 100

start_socs = [20, 20, 30, 20, 40, 20]
target_socs = [80, 70, 100, 90, 70, 70]
target_times = [16, 16, 15, 18, 20, 18]
bat_energies = [50, 50, 50, 50, 50, 50]
bus_volts = [400-i/2 for i in bus_lst]

# BEVs erzeugen
bev_list = []
for bus in bus_lst:
    bev = BEV(soc_start=start_socs[bus], soc_target=target_socs[bus],
              t_target=target_times[bus], e_bat=bat_energies[bus],
              resolution=resolution, bus_voltage=bus_volts[bus],
              home_bus=bus)
    bev_list.append(bev)

test = GLO(number_buses=buses, bevs=bev_list, resolution=resolution, s_trafo_kVA=s_trafo,
           households=None)

test.optimization_model.SOC.pprint()

# optimieren lassen
test.run_optimization_single_timestep(tee=True)

# Ergebnisse darstellen
test.plot_results()
