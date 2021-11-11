"""
Klasse GridLineOptimizer, die ein Model zur Optimierung der Ladeleistungen von Ladesäulen entlang eines Netzstrahls
erzeug. Um die Ergebnisse hinterher validieren zu können, ist auch direkt ein pandapower-Netz mit enthalten, mit dem
man nach der Optimierung die Ergebnisse überprüfen kann.
"""

_pandapower_available = True
_networkx_available = True

import pyomo.environ as pe

try:
    import pandapower as pp
    #raise ModuleNotFoundError

except ModuleNotFoundError:
    print('\nWARNING: module pandapower not available, some features are',
          'only available with pandapower\n')
    _pandapower_available = False

import matplotlib.pyplot as plt

try:
    import networkx as nx

except ModuleNotFoundError:
    print('\nWARNING: module networkx not available, some features are',
          'only available with networkx\n')
    _networkx_available = False


import time
from battery_electric_vehicle import BatteryElectricVehicle as BEV



class GridLineOptimizer:
    global _pandapwer_available
    global _networkx_available

    def __init__(self, number_buses, bev_buses, charger_locs=None, voltages=None, impedances=None,
                 resolution=60, s_trafo_kVA=100, solver='glpk'):
        self.current_timestep = 0
        self.number_buses = number_buses
        self.buses = self._make_buses()
        self.lines = self._make_lines()
        self.times = self._make_times()
        #self.buses_at_times = self._make_buses_at_times()
        if voltages == None:
            self.voltages = self._make_voltages()
        else:
            self.voltages = voltages
        self.i_max = 160
        self.u_min = 0.9*400/3**0.5
        self.s_trafo = s_trafo_kVA
        self.solver = solver
        self.solver_factory = pe.SolverFactory(self.solver)
        self.resolution = resolution
        if charger_locs == None:
            self.charger_locs = self.buses  # vielleicht lieber als dict: {1: True, 2: False...}?
        else:
            self.charger_locs = charger_locs

        if impedances == None:
            self.impedances = self._make_impedances()
        else:
            self.impedances = impedances

        self.bev_buses = bev_buses
        self.bevs = []
        self._make_bevs()

        self.optimization_model = self._setup_model()
        self.grid = self._setup_grid()

        # hier kommen dann die Ergebnisse für jeden Knoten zu jedem
        # timstep der Strom rein (die SOCs werden im BEV gespeichert)
        self.results_I = {bus: [] for bus in self.buses}


    def _make_buses(self):
        return list(range(self.number_buses))


    def _make_lines(self):
        return list(range(self.number_buses))

    # TODO: dafür sorgen, dass die Spanne von times gemäß der Auflösung und des horizonts angepasst wird
    def _make_times(self):
        return list(range(self.current_timestep, self.current_timestep+24))


    def _make_voltages(self):
        return {i: 400-i/2 for i in self.buses}


    #def _make_buses_at_times(self):
        #return {i: self.buses for i in self.times}


    def _make_impedances(self):
        return {i: 0.04 for i in self.lines}


    def _make_bevs(self):
        for bus in self.bev_buses:
            bev_bus_voltage = list(self.voltages)[bus]
            bev = BEV(home_bus=bus, e_bat=50, bus_voltage=bev_bus_voltage, resolution=self.resolution)
            print('BEV erzeugt an Bus', bus)
            self.bevs.append(bev)


    def _setup_model(self):
        # Model erzeugen
        model = pe.ConcreteModel('GridLineOptimization')

        # Sets als Indizes erzeugen
        model.buses = pe.Set(initialize=self.buses)
        model.lines = pe.Set(initialize=self.lines)
        model.times = pe.Set(initialize=self.times)

        # Parameter erzeugen
        model.impedances = pe.Param(model.lines, initialize=self.impedances)
        model.voltages = pe.Param(model.buses, initialize=self.voltages)
        model.u_min = self.u_min
        model.i_max = self.i_max

        # Entscheidungsvariablen erzeugen
        model.I = pe.Var(model.times*model.buses, domain=pe.PositiveReals, bounds=(0, 27))
        model.SOC = pe.Var(model.times*model.buses, domain=pe.PositiveReals, bounds=(0, 100))

        # Zielfunktion erzeugen
        def max_power_rule(model):
            return sum(sum(model.voltages[i]*model.I[j, i] for i in model.buses) for j in model.times)

        model.max_power = pe.Objective(rule=max_power_rule, sense=pe.maximize)


        # Einschränkungen festlegen
        def min_voltage_rule(model, t):
            return model.voltages[0] - sum(model.impedances[i] * sum(model.I[t, j] for j in range(i, len(model.buses)))
                                           for i in model.lines) >= model.u_min


        def max_current_rule(model, t):
            return sum(model.I[t, n] for n in model.buses) <= model.i_max


        def track_socs_rule(model, t, b):
            if t < self.current_timestep + 23:
                return (model.SOC[t, b] + model.I[t, b] * model.voltages[b] * self.resolution/60 / 1000
                        - model.SOC[t+1, b]) == 0

            else:
                return pe.Constraint.Skip


        model.min_voltage = pe.Constraint(model.times, rule=min_voltage_rule)
        model.max_current = pe.Constraint(model.times, rule=max_current_rule)
        model.track_socs = pe.Constraint(model.times*model.buses, rule=track_socs_rule)

        return model


    def _setup_grid(self):
        if not _pandapower_available:
            print('\nWARNING: unable to create grid\n')
            return None

        else:
            grid = pp.create_empty_network(name='OptimizationGrid')

            # Busse erzeugen
            pp.create_bus(grid, name='transformer mv', vn_kv=20)
            pp.create_bus(grid, name='transformer lv', vn_kv=0.4)

            for nr in range(self.number_buses):
                pp.create_bus(grid, name='bus '+str(nr), vn_kv=0.4)

            # Slack erzeugen
            pp.create_ext_grid(grid, bus=0)

            # Generator erzeugen
            pp.create_transformer_from_parameters(grid, hv_bus=0, lv_bus=1, sn_mva=self.s_trafo/1000,
                                                  vn_hv_kv=20, vn_lv_kv=0.4, vkr_percent=1.5, pfe_kw=0.4,
                                                  i0_percent=0.4, vk_percent=6)

            # Leitungen erzeugen
            for nr in range(self.number_buses-1):
                pp.create_line_from_parameters(grid, from_bus=nr+1, to_bus=nr+2, r_ohm_per_km=1,
                                               length_km=1/self.impedances[nr], name='line '+str(nr),
                                               x_ohm_per_km=0, c_nf_per_km=0, max_i_ka=0.142)

            return grid


    def display_target_function(self):
        self.optimization_model.max_power.pprint()


    def display_min_voltage_constraint(self):
        self.optimization_model.min_voltage.pprint()


    def display_max_current_constraint(self):
        self.optimization_model.max_current.pprint()


    def display_track_socs_constraint(self):
        self.optimization_model.track_socs.pprint()

    # nach jedem Optimierungsdurchlauf die Ergebnisse aus dem Model und
    # die SOCs in den BEVs speichern, I hier irgendwo...
    def _store_results(self):
        """
        Fragt die Optimierungsergebnisse der Entscheidungsvariablen I und SOC aus dem
        Optimierungsmodel ab und speichert diese.
        :return:
        """
        for num, bev in enumerate(self.bevs):
            # immer vom ersten (also aktuellen) timestep den entsprechenden SOC
            # wählen
            SOC = self.optimization_model.SOC[self.current_timestep, num].value
            bev.enter_soc(SOC)

        for bus in self.buses:
            # immer vom ersten (also aktuellen) timestep den entsprechenden I
            # wählen
            I = self.optimization_model.I[self.current_timestep, bus].value
            self.results_I[bus].append(I)


    # jetzt wo es run_optimization_rolling_horizon gibt, wird diese methode
    # eigentlich nicht mehr benötigt
    def run_optimization_single_timestep(self, **kwargs):
        self.solver_factory.solve(self.optimization_model, tee=kwargs['tee'])
        #return list(self.optimization_model.I[:, :].value)


    def run_optimization_rolling_horizon(self, complete_horizon, **kwargs):
        for i in range(complete_horizon):
            print('aktueller Zeitschritt:', i)
            # Modell entsprechend neu aufbauen (times geht von i bis i+24)
            self.current_timestep = i
            self.times = self._make_times()
            self.optimization_model = self._setup_model()
            self.solver_factory.solve(self.optimization_model, tee=kwargs['tee'])
            self._store_results()


    def plot_grid(self):
        if not _networkx_available:
            print('\nWARNING: unable to plot grid\n')

        else:
            # simple_plotly und simple_plot gehen nicht
            graph = nx.DiGraph()
            edges = set([(i, i+1) for i in range(len(self.buses)-1)])
            graph.add_nodes_from(list(self.buses))
            graph.add_edges_from(list(edges))
            pos = ({i: (i, 4) for i in self.buses})

            fig, ax = plt.subplots(figsize=(12, 12))

            nx.draw_networkx_nodes(graph, pos=pos, ax=ax, node_color='lightgray',
                                   edgecolors='black', node_size=2000)

            #nx.draw_networkx_labels(graph, pos=pos, ax=ax, labels=dict(zip(nodes, nodes)),
                                    #font_size=20)

            nx.draw_networkx_edges(graph, pos=pos, ax=ax, node_size=2000, arrowsize=25)

            #nx.draw_networkx_edge_labels(graph, pos=pos, ax=ax, edge_labels=distances,
                                         #font_size=16, rotate=False)

            plt.show()


    def list_bevs(self):
        for bev in self.bevs:
            print('---------------------------------')
            print(f'BEV an Bus {bev.home_bus}')
            print(f'mit Batterie {bev.e_bat} kWh')
            #print(f'und SOC {bev.current_soc} %')
            # TODO noch dafür sorgen, dass jedes BEV die korrekte Knotenspannung bekommt
            print(f'an Knotenspannung {bev.bus_voltage} V')


    def plot_results(self):
        pass




if __name__ == '__main__':
    t0 = time.time()
    test = GridLineOptimizer(6, bev_buses=list(range(6)), resolution=60)
    print(test.bevs)
    test.list_bevs()


    # print(test.buses)
    # print(test.lines)
    # print(test.impedances)
    # print(test.u_min)
    test.display_target_function()
    test.display_min_voltage_constraint()
    test.display_max_current_constraint()
    test.display_track_socs_constraint()
    #res = test.run_optimization_single_timestep(tee=True)
    # #for i, val in enumerate(res):
    #     #print(f'Strom am Knoten {i}: {val}')
    # #
    #dt = time.time() - t0
    #print('Laufzeit', dt)
    print('starte rolling horizon: \n')
    test.run_optimization_rolling_horizon(24, tee=False)
    dt = time.time() - t0
    print('fertig nach', dt, 'sek')
    # test.optimization_model.I.pprint()
    # print(res)
    #test.optimization_model.SOC.pprint()
    test.optimization_model.I.pprint()
    print(test.results_I[0])


