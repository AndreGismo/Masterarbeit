"""
Haushalts-Klasse in der die Informationen zu den Lastprofilen der Haushalte steckt
"""
import pandas as pd
import matplotlib.pyplot as plt

class Household:
    _data_source = '../Data/Nuernberg_SLP.csv'
    _e_norm = 3000  # kWh jährlicher Verbrauch
    _res_norm = '15min'  # min Auflösung
    _len_norm = 35040

    def __init__(self, home_bus, annual_demand=3000, resolution=15):
        self.home_bus = home_bus
        self.annual_demand = annual_demand
        self.resolution = str(resolution)+'min'
        self.load_profile = None
        self.calc_load_profile()


    def calc_load_profile(self):
        self.load_profile = pd.read_csv(type(self)._data_source)
        self.load_profile *= self.annual_demand/type(self)._e_norm
        self.load_profile.index = pd.date_range('2021', periods=type(self)._len_norm, freq=type(self)._res_norm)
        if not self.resolution == type(self)._res_norm:
            #wenn die gewünschte Auflösung eine andere ist als die gegebene, dann resamplen
            self.load_profile = self.load_profile.resample(self.resolution).mean()


    def plot_load_profile(self):
        self.load_profile.plot()
        plt.show()



if __name__ == '__main__':
    slp = Household(5, 5000, resolution=60)
    slp.plot_load_profile()




