"""
Klasse zur Simulation von Battery Electric Vehicles (BEVs). Der SOC wird getrackt und dann an die Optimierung
zurückgemeldet. Außerdem wird festgehalten, an welcher Ladesäule (also an welchem Bus im Netz) das BEV lädt und
wie groß die Batterie ist.
"""

class BatteryElectricVehicle:
    def __init__(self, home_bus, e_bat=50, soc_start=50, soc_target=100, t_target=17,
                 t_start=14, resolution=60, current_timestep=0, recurring='daily'):
        self.home_bus = home_bus
        self.e_bat = e_bat
        #self.bus_voltage = bus_voltage
        self.soc_start = soc_start
        self.soc_target = soc_target
        self.resolution = resolution
        self.t_target = int(t_target * 60 / self.resolution)
        self.t_start = int(t_start * 60/self.resolution)
        self.soc_list = [soc_start]
        self.is_loading = True
        self.current_timestep = current_timestep
        self.recurring = recurring
        self.horizon_width = None # bekommt von GLO mitgeteilt
        self.occupancies = None # wird auch von GLO aus aufgerufen


    # wird von GLO aus aufgerufen
    def make_occupancies(self):
        occupancies = [False for _ in range(int(self.horizon_width * 60/self.resolution))]
        offset = 24 * 60/self.resolution
        for num in range(int(self.horizon_width/24)):
            occupancies[int(self.t_start+num*offset):int(self.t_target+num*offset)] = [True for _ in range(self.t_target-self.t_start)]
        self.occupancies = occupancies


    # wird von GLO aus aufgerufen
    def set_horizon_width(self, width_hrs):
        self.horizon_width = width_hrs


    def enter_soc(self, soc):
        self.soc_list.append(soc)
        # wenn der Ladestand größer soc_target, dann "fährt der weg"
        if soc >= self.soc_target:
            self.is_loading = False


    def plot_soc(self):
        pass


    def check_availability(self, timestep):
        self.current_timestep = timestep
        if self.current_timestep < self.t_start:
            self.is_loading = False
            return False
        elif self.current_timestep >= self.t_start and self.current_timestep <= self.t_target:
            self.is_loading = True
            return True
        else:
            self.is_loading = False
            return False
