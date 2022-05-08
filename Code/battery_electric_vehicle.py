"""
Klasse zur Simulation von Battery Electric Vehicles (BEVs). Der SOC wird getrackt und dann an die Optimierung
zurückgemeldet. Außerdem wird festgehalten, an welcher Ladesäule (also an welchem Bus im Netz) das BEV lädt und
wie groß die Batterie ist.
"""
import numpy as np

class BatteryElectricVehicle:
    def __init__(self, home_bus, e_bat=50, soc_start=50, soc_target=100, t_target=17,
                 t_start=14, resolution=60, current_timestep=0, p_load=11,
                 recurring='daily'):
        self.home_bus = home_bus
        self.e_bat = e_bat
        #self.bus_voltage = bus_voltage
        self.soc_start = soc_start
        self.soc_target = soc_target
        self.resolution = resolution
        self.t_target = int(t_target * 60 / self.resolution)
        self.t_start = int(t_start * 60/self.resolution)
        self.p_load = p_load
        self.soc_list = [soc_start]
        self.is_loading = True
        self.current_timestep = current_timestep
        self.recurring = recurring
        self.horizon_width = None # bekommt von GLO mitgeteilt
        self.occupancies = None # wird auch von GLO aus aufgerufen
        self.current_soc = soc_start
        self._make_waiting_times()


    def update_soc(self, value):
        self.current_soc = value
        #print(f'SOC of BEV at node{self.home_bus} at timestep {self.current_timestep}: {self.current_soc} %')


    def get_current_power(self, timestep, cap):
        if timestep < self.t_start:
            return 0
        elif timestep >= self.t_start  and timestep < self.t_target:
            if self.current_soc <= 80:
                self.calc_new_soc(self.p_load*cap)
                return self.p_load*cap
            elif self.current_soc > 80 and self.current_soc <= 100:
                # calculate according to exponential decrease formula
                p_load_calc = self.calc_p_load()
                self.calc_new_soc(p_load_calc*cap)
                return p_load_calc*cap
            else:
                return 0
        else:
            return 0


    def calc_new_soc(self, power):
        self.current_soc += (power * self.resolution/60)/self.e_bat*100
        if self.current_soc > 100:
            self.current_soc = 100
        print(f'current soc of BEV at node {self.home_bus}: {self.current_soc}')


    def calc_p_load(self):
        # calculate load stop power
        p_ls = 4.2/3.9*0.03*self.e_bat
        # calculate loading correcture factor
        kl = (100-80)/np.log(self.p_load/p_ls)
        # calculate P(SOC)
        p_soc = self.p_load * np.exp((80-self.current_soc)/kl)
        return p_soc


    def _make_waiting_times(self):
        self.waiting_times = [0 for _ in range(int(24*60/self.resolution))]
        waited = 0
        for num, val in enumerate(self.waiting_times):
            if num >= self.t_start and num < self.t_target:
                self.waiting_times[num] = waited
                waited += 1


    # wird von GLO aus aufgerufen
    # deprecated
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

    # deprecated
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
