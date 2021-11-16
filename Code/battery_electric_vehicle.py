"""
Klasse zur Simulation von Battery Electric Vehicles (BEVs). Der SOC wird getrackt und dann an die Optimierung
zurückgemeldet. Außerdem wird festgehalten, an welcher Ladesäule (also an welchem Bus im Netz) das BEV lädt und
wie groß die Batterie ist.
"""

class BatteryElectricVehicle:
    def __init__(self, home_bus, e_bat, bus_voltage, soc_start=50, soc_target=100, t_target=17, resolution=60):
        self.home_bus = home_bus
        self.e_bat = e_bat
        self.bus_voltage = bus_voltage
        self.soc_start = soc_start
        self.soc_target = soc_target
        self.resolution = resolution
        self.t_target = int(t_target * 60 / self.resolution)
        self.soc_list = [soc_start]
        self.is_loading = True


    def enter_soc(self, soc):
        self.soc_list.append(soc)
        # wenn der Ladestand größer soc_target, dann "fährt der weg"
        if soc >= self.soc_target:
            self.is_loading = False


    def plot_soc(self):
        pass
