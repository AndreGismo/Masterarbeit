"""
Author: André Ulrich
--------------------
Class for simulating BEVs.

Usage: Just create them and pass them to the constructor of GridLineOptimizer.
They can also be passed to Simulation_Handler.run_unoptimized_sim. In this case
they use the implemented P(SOC) characteristic to calculate the charging power.

Version history (only the most relevant points, full history is available on github):
-------------------------------------------------------------------------------------------------
V.1: first working description of a BEV

V.2: added functionality for charging according to P(SOC) characteristic

all the other commits in much more detail are available here:
https://github.com/AndreGismo/Masterarbeit/tree/submission)
"""
import numpy as np

class BatteryElectricVehicle:
    """
    class to simulate BEVs
    """
    def __init__(self, home_bus, soc_start, soc_target, t_target,
                 t_start, resolution, current_timestep=0, p_load=11,
                 e_bat=50, recurring='daily'):
        """
        Build BEV

        :param home_bus: bus of the grid line (including 0) where the BEV is charging
        :param soc_start: SOC at start of charging
        :param soc_target: desired SOC at end of charging
        :param t_target: deadline for charging to be finished
        :param t_start: charging start time
        :param resolution: resolution in time
        :param current_timestep: the current time step
        :param p_load: nominal loading power of the charging station at home_bus
        :param e_bat: nominal capacity of the BEVs battery
        :param recurring: wheather or not the BEV starts charging at multiple days (currently not really used)
        """
        self.home_bus = home_bus
        self.e_bat = e_bat
        self.soc_start = soc_start
        self.soc_target = soc_target
        self.resolution = resolution
        self.t_target = int(t_target * 60 / self.resolution)
        self.t_start = int(t_start * 60/self.resolution)
        self.p_load = p_load
        self.soc_list = []
        self.is_loading = True
        self.current_timestep = current_timestep
        self.recurring = recurring
        self.horizon_width = None
        self.current_soc = soc_start


    def update_soc(self, value):
        """
        updates the current_soc of the BEV with the value from the optimization.

        :param value: SOC [%]
        :return: None
        """
        self.current_soc = value
        self.soc_list.append(value)


    def get_current_power(self, timestep, cap):
        """
        calculate the loading power of the BEV at the current timestep.

        :param timestep: the current timestep
        :param cap: power cap of the P(U)-controller
        :return: load power of BEV including possible
        cap of the controller
        """
        if timestep < self.t_start:
            # no loading if BEV is not at charger
            return 0

        elif timestep >= self.t_start  and timestep < self.t_target:
            if self.current_soc <= 80:
                # load power = nominal power
                self.calc_new_soc(self.p_load*cap)
                return self.p_load*cap

            elif self.current_soc > 80 and self.current_soc <= 100:
                # calculate according to exponential decrease formula
                p_load_calc = self.calc_p_load()
                self.calc_new_soc(p_load_calc*cap)
                return p_load_calc*cap

            else:
                # no more loading if SOC = 100%
                return 0

        else:
            # no loading if BEV is not at charger
            return 0


    def calc_new_soc(self, power):
        """
        calculate new soc fo the timestep according to the load power
        and battery size.

        :param power: load power of the current timestep
        :return: None
        """
        self.current_soc += (power * self.resolution/60)/self.e_bat*100
        if self.current_soc > 100:
            self.current_soc = 100
        #print(f'current soc of BEV at node {self.home_bus}: {self.current_soc}')


    def calc_p_load(self):
        """
        calculate new load power for the timestep according to P(SOC)
        characteristic

        :return: new load power
        """
        # calculate load stop power
        p_ls = 4.2/3.9*0.03*self.e_bat
        # calculate loading correcture factor
        kl = (100-80)/np.log(self.p_load/p_ls)
        # calculate P(SOC)
        p_soc = self.p_load * np.exp((80-self.current_soc)/kl)
        return p_soc


    def reset_soc(self):
        """
        resets the BEVs SOC (only needed, in case a unoptimized grid simultion
        is run to make sure the calculation of charging power according to
        P(SOC) characteristic goes well (because a rolling horizon optimization
        beforhand might have already tempered the SOCs)).
        Only gets called from inside Simulation_Handler.run_unoptimized_sim.

        :return: None
        """
        self.current_soc = self.soc_start


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
