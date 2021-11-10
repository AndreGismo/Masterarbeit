"""
Klasse zur Simulation von Battery Electric Vehicles (BEVs). Der SOC wird getrackt und dann an die Optimierung
zurückgemeldet. Außerdem wird festgehalten, an welcher Ladesäule (also an welchem Bus im Netz) das BEV lädt und
wie groß die Batterie ist.
"""

class BatteryElectricVehicle:
    def __init__(self, home_bus, e_bat, bus_voltage, soc_start=50, t_occupancy=None, resolution=None):
        self.home_bus = home_bus
        self.e_bat = e_bat
        self.bus_voltage = bus_voltage
        self.soc_start = soc_start
        self.t_occupancy = t_occupancy
        self.resolution = resolution
        self.current_soc = soc_start


    def plot_soc(self):
        pass
