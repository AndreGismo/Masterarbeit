"""
Haushalts-Klasse in der die Informationen zu den Lastprofilen der Haushalte steckt
"""
import pandas as pd
import matplotlib.pyplot as plt

class Household:
    """
    class to simulate household load profiles
    """
    _data_source = '../Data/Nuernberg_SLP.csv'
    _e_norm = 3000  # kWh jährlicher Verbrauch
    _res_norm = '15min'  # min Auflösung
    _len_norm = 35040

    def __init__(self, home_bus, annual_demand=3000, resolution=15):
        """
        create household
        :param home_bus: bus of the grid line (including 0) where the household is attached
        :param annual_demand: annual energy demand (kWh) for scaling the load profile
        :param resolution: resolution in time (minutes) for the load profile
        """
        self.home_bus = home_bus
        self.annual_demand = annual_demand
        self.resolution = str(resolution)+'min'
        self.load_profile = None
        self.calc_load_profile()


    def calc_load_profile(self):
        self.load_profile = pd.read_csv(type(self)._data_source)
        # scale according to annual demand
        self.load_profile *= self.annual_demand/type(self)._e_norm
        # datetime index for easy calculation of mean/interpolate
        self.load_profile.index = pd.date_range('2021', periods=type(self)._len_norm, freq=type(self)._res_norm)
        asked_res = int(self.resolution.rstrip('min'))
        nom_res = int(type(self)._res_norm.rstrip('min'))

        # calculate mean if the asked res is greater than nominal resolution
        if asked_res >= nom_res:
            self.load_profile = self.load_profile.resample(self.resolution).mean()

        # calculate interpolated points if the asked res is smaller than nominal resolution
        else:
            self.load_profile = self.load_profile.resample(self.resolution).interpolate()

        self.load_profile = self.load_profile.values


    def plot_load_profile(self):
        plt.plot(list(range(len(self.load_profile))), self.load_profile)
        plt.show()


    def raise_demand(self, start, end, demand, recurring=None):
        """
        alter the load profile of the household
        :param start: start time of additional load
        :param end: end time of additional load
        :param demand: amount of additional load (kW), negative values lower the demand
        :param recurring: wheather the altered timerange should be repeated
        :return:
        """
        if recurring == 'daily':
            cycles = 364
            offset = int(24 * 60/int(self.resolution.rstrip('min'))) # 24 hrs/day
        elif recurring == 'weekly':
            cycles = 51
            offset = int(168 * 60/int(self.resolution.rstrip('min'))) # 168 hrs/week

        start = int(start * 60/int(self.resolution.rstrip('min')))
        end = int(end * 60/int(self.resolution.rstrip('min')))

        if recurring == None:
                self.load_profile[start:end] += demand

        else:
            for cycle in range(cycles):
                # change at multiple time intervals
                self.load_profile[int(start+cycle*offset):int(end+cycle*offset)] += demand



if __name__ == '__main__':
    slp = Household(5, 5000, resolution=5)
    #slp.raise_demand(10, 14, 1200, recurring='weekly')
    slp.plot_load_profile()




