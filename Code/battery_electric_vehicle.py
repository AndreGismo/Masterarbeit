"""
Klasse zur Simulation von Battery Electric Vehicles (BEVs). Der SOC wird getrackt und dann an die Optimierung
zurückgemeldet. Außerdem wird festgehalten, an welcher Ladesäule (also an welchem Bus im Netz) das BEV lädt und
wie groß die Batterie ist.
"""

class BatteryElectricVehicle:
    def __init__(self, home_node, e_bat, t_occupancy):
        self.home_node = home_node
        self.e_bat = e_bat
        self.t_occupancy = t_occupancy