"""
Klasse GridLineOptimizer, die ein Model zur Optimierung der Ladeleistungen von Ladesäulen entlang eines Netzstrahls
erzeug. Um die Ergebnisse hinterher validieren zu können, ist auch direkt ein pandapower-Netz mit enthalten, mit dem
man nach der Optimierung die Ergebnisse überprüfen kann.
"""

import pyomo.environ as pe
import pandapower as pp
import matplotlib.pyplot as plt
import networkx as nx
import time



class GridLineOptimizer:
    def __init__(self, number_buses, charger_locs=None, voltages=None, impedances=None,
                 s_trafo_kVA=100, solver='glpk'):
        self.number_buses = number_buses
        self.buses = self._make_buses()
        self.lines = self._make_lines()
        self.times = set(range(24))
        self.buses_at_times = self._make_buses_at_times()
        if voltages == None:
            self.voltages = self._make_voltages()
        else:
            self.voltages = voltages
        self.i_max = 160
        self.u_min = 0.9*400/3**0.5
        self.s_trafo = s_trafo_kVA
        self.solver = solver
        self.solver_factory = pe.SolverFactory(self.solver)
        if charger_locs == None:
            self.charger_locs = self.buses  # vielleicht lieber als dict: {1: True, 2: False...}?
        else:
            self.charger_locs = charger_locs

        if impedances == None:
            self.impedances = self._make_impedances()
        else:
            self.impedances = impedances

        self.optimization_model = self._setup_model()
        self.grid = self._setup_grid()


    def _make_buses(self):
        return set(range(self.number_buses))


    def _make_lines(self):
        return set(range(self.number_buses))


    def _make_voltages(self):
        return {i: 400 for i in self.buses}


    def _make_buses_at_times(self):
        return {i: self.buses for i in self.times}


    def _make_impedances(self):
        return {i: 0.04 for i in self.lines}


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
        model.I.pprint()

        # Zielfunktion erzeugen
        def max_power_rule(model):
            return sum(sum(model.voltages[i]*model.I[j, i] for i in model.buses) for j in model.times)

        model.max_power = pe.Objective(rule=max_power_rule, sense=pe.maximize)

        # TODO: das verursacht den Fehler beim Constraint (weil I jetzt natürlich falsch indexiert ist => anpassen!
        # Einschränkungen festlegen
        def min_voltage_rule(model, t):
            return model.voltages[0] - sum(model.impedances[i] * sum(model.I[t, j] for j in range(i, len(model.buses)))
                                           for i in model.lines) >= model.u_min


        def max_current_rule(model, t):
            return sum(model.I[t, n] for n in model.buses) <= model.i_max


        model.min_voltage = pe.Constraint(model.times, rule=min_voltage_rule)
        model.max_current = pe.Constraint(model.times, rule=max_current_rule)

        return model


    def _setup_grid(self):
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


    def run_optimization_single_timestep(self, **kwargs):
        self.solver_factory.solve(self.optimization_model, tee=kwargs['tee'])
        return list(self.optimization_model.I[:, :].value)


    def plot_grid(self):
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


if __name__ == '__main__':
    t0 = time.time()
    test = GridLineOptimizer(8)


    print(test.buses)
    print(test.lines)
    print(test.impedances)
    print(test.u_min)
    test.display_target_function()
    test.display_min_voltage_constraint()
    test.display_max_current_constraint()
    res = test.run_optimization_single_timestep(tee=True)
    #for i, val in enumerate(res):
        #print(f'Strom am Knoten {i}: {val}')
    #
    dt = time.time() - t0
    print('Laufzeit', dt)
    #test.plot_grid()
    test.optimization_model.I.pprint()
    print(res)