"""
Klasse GridLineOptimizer, die ein Model zur Optimierung der Ladeleistungen von Ladesäulen entlang eines Netzstrahls
erzeug. Um die Ergebnisse hinterher validieren zu können, ist auch direkt ein pandapower-Netz mit enthalten, mit dem
man nach der Optimierung die Ergebnisse überprüfen kann.

Eventuell macht die rolling horizon Betrachtung hier auch gar keinen Sinn? Weil sich ja noch keine Werte im Lauf der
Zeit mal zufällig ändern (wie das in der Realität der Fall wäre). Vielleicht genügt ja ein Durchlauf mit dem 24std
Horizont - oder der Fehler liegt einfach in der Logik von store_results

Vielleicht Weite des horizons auch als parameter mit aufnehmen

Beim upper und lower bounds von SOC und I
"""

_pandapower_available = True
_networkx_available = True
_pandas_available = True
_matplotlib_available = True

import pyomo.environ as pe

try:
    import pandas as pd

except ModuleNotFoundError:
    print('\nWARNING: module pandas not available, some features are',
          'only available with pandas\n')
    _pandas_available = False

try:
    import pandapower as pp
    #raise ModuleNotFoundError

except ModuleNotFoundError:
    print('\nWARNING: module pandapower not available, some features are',
          'only available with pandapower\n')
    _pandapower_available = False

try:
    import matplotlib.pyplot as plt

except ModuleNotFoundError:
    print('\nWARNING: module matplotlib not available, some features are',
          'only available with matplotlib\n')
    _matplotlib_available = False

try:
    import networkx as nx

except ModuleNotFoundError:
    print('\nWARNING: module networkx not available, some features are',
          'only available with networkx\n')
    _networkx_available = False

import time
from battery_electric_vehicle import BatteryElectricVehicle as BEV
from household import Household as HH



class GridLineOptimizer:
    global _pandapwer_available
    global _networkx_available
    global _pandas_available
    global _matplotlib_available

    def __init__(self, number_buses, bevs, households, charger_locs=None, horizon_width=24,
                 voltages=None, impedances=None, resolution=60, s_trafo_kVA=100, solver='glpk'):
        self.current_timestep = 0
        self.resolution = resolution
        self.horizon_width = horizon_width
        self.number_buses = number_buses
        self.buses = self._make_buses()
        self.lines = self._make_lines()
        #self.times = self._make_times()
        self._make_times()
        if voltages == None:
            self.voltages = self._make_voltages()
        else:
            self.voltages = voltages
        #self.i_max = 160   # 160
        self.u_min = 0.9*400
        self.s_trafo = s_trafo_kVA
        self.i_max = s_trafo_kVA*1000 / 400
        self.solver = solver
        self.solver_factory = pe.SolverFactory(self.solver)
        #self.charger_locs = None

        if impedances == None:
            self.impedances = self._make_impedances()
        else:
            self.impedances = impedances

        self._make_bev_dict(bevs)
        self._setup_bevs()
        self.households = households
        #self._determine_charger_locs()

        self._make_soc_lower_bounds()
        self._make_soc_upper_bounds()

        self.prepare_i_lower_bounds()
        self.prepare_i_upper_bounds()

        self.optimization_model = self._setup_model()
        self.grid = self._setup_grid()

        # hier kommen dann die Ergebnisse für jeden Knoten zu jedem
        # timstep der Strom rein (die SOCs werden im BEV gespeichert)
        self.results_I = {bus: [] for bus in self.buses}


    # wird nur bei dem allerersten Durchlauf der Optimierung genutzt
    def _make_soc_lower_bounds(self):
        soc_lower_bounds = {bev.home_bus: [bev.soc_start for _ in range(len(self.times))] for bev in self.bevs.values()}
        for bev in self.bevs.values():
            # dafür sorgen, dass an demjenigen Zeitpunkt, wo die geladen sein wollen t_target
            # der gewünschte Ladestand soc_target dasteht
            soc_lower_bounds[bev.home_bus][bev.t_target] = bev.soc_target
        self.soc_lower_bounds = soc_lower_bounds


    # wird nur bei dem allerersten Durchlauf der Optimierung genutzt
    def _make_soc_upper_bounds(self):
        soc_upper_bounds = {bev.home_bus: [bev.soc_target for _ in range(len(self.times))] for bev in self.bevs.values()}
        for bev in self.bevs.values():
            # dafür sorgen, dass beim Startzeitpunkt die upper bound gleich der lower bound
            # (also soc start) ist (bei anderen Startpunkten als 0 noch entsprechendes
            # t_start in BEV einführen und hier statt 0 nutzen
            soc_upper_bounds[bev.home_bus][0] = bev.soc_start
        self.soc_upper_bounds = soc_upper_bounds


    def prepare_i_lower_bounds(self):
        i_lower_bounds = {bev.home_bus: [0 for _ in range(len(self.times))] for bev in self.bevs.values()}
        self.i_lower_bounds = i_lower_bounds


    def prepare_i_upper_bounds(self):
        i_upper_bounds = {bev.home_bus: [27 for _ in range(len(self.times))] for bev in self.bevs.values()}
        for bev in self.bevs.values():
            i_upper_bounds[bev.home_bus][0:bev.t_start] = [0 for _ in range(bev.t_start)]
        self.i_upper_bounds = i_upper_bounds


    def _make_bev_dict(self, bevs):
        bev_dict = {bev.home_bus: bev for bev in bevs}
        self.bevs = bev_dict


    def _setup_bevs(self):
        for bev in self.bevs.values():
            bev.set_horizon_width(self.horizon_width)
            bev.make_occupancies()


    def _prepare_soc_lower_bounds(self):
        soc_lower_bounds = {bev.home_bus: [bev.soc_start for _ in range(len(self.times))] for bev in self.bevs.values()}
        self.soc_lower_bounds = soc_lower_bounds


    def _prepare_soc_upper_bounds(self):
        soc_upper_bounds = {bev.home_bus: [bev.soc_target for _ in range(len(self.times))] for bev in self.bevs.values()}
        self.soc_upper_bounds = soc_upper_bounds


    def _make_buses(self):
        return list(range(self.number_buses))


    def _make_lines(self):
        return list(range(self.number_buses))


    def _make_times(self):
        self.times = list(range(self.current_timestep, self.current_timestep+self.horizon_width*int(60/self.resolution)))


    def _make_voltages(self):
        return {i: 400-i/2 for i in self.buses}


    def _make_impedances(self):
        return {i: 0.04 for i in self.lines}  # 0.04


    # wird nicht mehr benutzt, da BEVs außerhalb erzeugt werden
    def _make_bevs(self):
        for bus in self.bev_buses:
            bev_bus_voltage = self.voltages[bus]
            start_soc = bus * 5 + 5
            bev = BEV(home_bus=bus, e_bat=50, bus_voltage=bev_bus_voltage, resolution=self.resolution,
                      soc_start=start_soc)
            self.bevs.append(bev)


    # wird nicht mehr benötigt, da das jetzt direkt über self.bevs ersichtlich ist
    def _determine_charger_locs(self):
        locations = []
        for bev in self.bevs:
            locations.append(bev.home_bus)
        self.charger_locs = locations


    def _prepare_next_timestep(self):
        self.current_timestep += 1
        self._make_times()
        # Ergebnisse des zweiten timesteps des vorherigen horizons nehmen und an
        # die erste Stelle der soc_upper und lower_bounds schreiben (dasselbe auch
        # für I?)
        SOCs2 = []
        for bus in self.bevs:
            SOCs2.append(self.optimization_model.SOC[self.current_timestep, bus].value)

        self._prepare_soc_lower_bounds()
        self._prepare_soc_upper_bounds()

        for num, bev in enumerate(self.bevs.values()):
            # an der ersten Stelle als lower und upper bound jetzt das Ergebnis
            # schreiben
            self.soc_lower_bounds[bev.home_bus][0] = SOCs2[num]
            self.soc_upper_bounds[bev.home_bus][0] = SOCs2[num]

        # Ergebnisse des zweiten timesteps des vorherigen horizons nehmen und für
        # i upper und lower_bounds an die erste Stelle schreiben
        Is2 = []
        for bus in self.bevs:
            Is2.append(self.optimization_model.I[self.current_timestep, bus].value)

        self.prepare_i_lower_bounds()
        self.prepare_i_upper_bounds()

        for num, bev in enumerate(self.bevs.values()):
            self.i_lower_bounds[bev.home_bus][0] = Is2[num]
            self.i_upper_bounds[bev.home_bus][0] = Is2[num]





    def _setup_model(self):
        # Model erzeugen
        model = pe.ConcreteModel('GridLineOptimization')

        # Sets als Indizes erzeugen
        model.buses = pe.Set(initialize=self.buses)
        model.lines = pe.Set(initialize=self.lines)
        model.times = pe.Set(initialize=self.times)
        model.charger_buses = pe.Set(initialize=[bev.home_bus for bev in self.bevs.values()])

        # Parameter erzeugen
        model.impedances = pe.Param(model.lines, initialize=self.impedances)
        model.voltages = pe.Param(model.buses, initialize=self.voltages)
        model.u_min = self.u_min
        model.i_max = self.i_max

        def get_household_currents(model, time, bus):
            # getielt durch die Spannung an dem Knoten, weil es ja Strom sein soll
            return self.households[bus].load_profile[time] / self.voltages[bus]

        model.household_currents = pe.Param(model.times*model.buses, initialize=get_household_currents)

        # Entscheidungsvariablen erzeugen (dafür erstmal am besten ein array (timesteps x buses)
        # wo überall nur 50 drinsteht (oder was man dem BEV halt als coc_start übergeben hatte))
        # erzeugen und diese Werte als lower bound ausgeben
        def get_soc_bounds(model, time, bus):
            return (self.soc_lower_bounds[bus][time-self.current_timestep], self.soc_upper_bounds[bus][time-self.current_timestep])

        def get_i_bounds(model, time, bus):
            return (self.i_lower_bounds[bus][time-self.current_timestep], self.i_upper_bounds[bus][time-self.current_timestep])

        model.I = pe.Var(model.times*model.charger_buses, domain=pe.PositiveReals, bounds=get_i_bounds)
        model.SOC = pe.Var(model.times*model.charger_buses, domain=pe.PositiveReals, bounds=get_soc_bounds)

        # Zielfunktion erzeugen
        def max_power_rule(model):
            #return sum(sum(model.voltages[i]*model.I[j, i] for i in model.charger_buses) for j in model.times)
            return sum(sum(model.SOC[t+1, b] - model.SOC[t, b] for b in model.charger_buses if t < len(model.times)-1) for t in model.times)
            #return sum(sum(model.SOC[t, b] - model.SOC[model.times.prevw(t), b] for b in model.charger_buses) for t in model.times)


        model.max_power = pe.Objective(rule=max_power_rule, sense=pe.maximize)

        # Einschränkungen festlegen
        def min_voltage_rule(model, t):
            return model.voltages[0] - sum(model.impedances[i] * (sum(model.household_currents[t, j] for j in model.buses if j > i)
                                                                  +sum(model.I[t, j] for j in model.charger_buses if j > i))
                                           for i in model.lines) >= model.u_min


        def max_current_rule(model, t):
            return sum(model.I[t, b] for b in model.charger_buses) + sum(model.household_currents[t, b] for b in model.buses) <= model.i_max


        def track_socs_rule(model, t, b):
            # schauen, dass man immer nur bis zum vorletzten timestep geht (weil es
            # sonst kein t+1 mehr geben würde beim letzten timestep)
            if t < self.current_timestep + self.horizon_width*60/self.resolution-1:#23:
                return (model.SOC[t, b] + model.I[t, b] * model.voltages[b] * self.resolution/60 / 1000
                        / self.bevs[b].e_bat*100 - model.SOC[t+1, b]) == 0

            else:
                return pe.Constraint.Skip

        #TODO: schauen, dass wenn der rolling horizon so weit fortgeschritten ist, dass der Ziel-
        # Zeitpunkt nicht mehr im aktuell betrachteten Horizont enthalten ist, dass dann dieser
        # Constraint auch nicht mehr auftaucht (geht vielleicht schon automatisch durch t == t_target)
        def ensure_final_soc_rule(model,  b):
            t_end = self.bevs[b].t_target
            if t_end - self.current_timestep > 0:
                return sum(model.voltages[b] * model.I[t, b] for t in range(self.current_timestep, t_end))* self.resolution/60\
                /1000 / self.bevs[b].e_bat * 100 <= (self.bevs[b].soc_target - self.bevs[b].soc_list[self.current_timestep-1])#- self.bevs[b].soc_start)
            else:
                return pe.Constraint.Skip


        model.min_voltage = pe.Constraint(model.times, rule=min_voltage_rule)
        model.max_current = pe.Constraint(model.times, rule=max_current_rule)
        model.track_socs = pe.Constraint(model.times*model.charger_buses, rule=track_socs_rule)
        # mit diesem Constraint kommt dasselbe raus, als hätte man nur track_socs aktiv
        #model.ensure_final_soc = pe.Constraint(model.buses, rule=ensure_final_soc_rule)

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


    def display_ensure_final_soc_constraint(self):
        self.optimization_model.ensure_final_soc.pprint()


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


    # simuliert einen Tag mit den Werten für einen Tag, fixer Horizont
    def run_optimization_single_timestep(self, **kwargs):
        self.solver_factory.solve(self.optimization_model)#, tee=kwargs['tee'])
        #return list(self.optimization_model.I[:, :].value)


    def run_optimization_rolling_horizon(self, complete_horizon, **kwargs):
        steps = int(complete_horizon * 60/self.resolution)
        for _ in range(steps):
            print(self.times)
            self.optimization_model.max_power.pprint()
            self.run_optimization_single_timestep()
            self._prepare_next_timestep()
            self.optimization_model = self._setup_model()


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
            print(f'an Knotenspannung {bev.bus_voltage} V')


    def plot_results(self):
        if not _pandas_available or not _matplotlib_available:
            print('\nWARNING: unable to plot results\n')

        else:
            # erstmal die ergebnisse aus dem Modell abfragen
            SOCs = {bus: [] for bus in self.bevs}
            for time in self.times:
                for bus in self.bevs:
                    SOCs[bus].append(self.optimization_model.SOC[time, bus].value)

            Is = {bus: [] for bus in self.bevs}
            for time in self.times:
                for bus in self.bevs:
                    Is[bus].append(self.optimization_model.I[time, bus].value)

            SOCs_df = pd.DataFrame(SOCs)
            SOCs_df.index = pd.date_range(start='2021', periods=len(SOCs_df), freq=str(self.resolution)+'min')
            Is_df = pd.DataFrame(Is)
            Is_df.index = pd.date_range(start='2021', periods=len(SOCs_df), freq=str(self.resolution) + 'min')

            fig, ax = plt.subplots(2, 1, figsize=(15, 15), sharex=False)
            for column in SOCs_df.columns:
                ax[0].plot(SOCs_df.index, SOCs_df[column], marker='o', label=f'SOC der Batterie am Knoten {column}')
            ax[0].legend()
            ax[0].grid()
            ax[0].set_ylabel('SOC [%]')
            ax[0].set_title('SOC über der Zeit - Ergebnisse der Optimierung')

            for column in Is_df.columns:
                ax[1].plot(Is_df.index, Is_df[column], marker='o', label=f'Strom in die Batterie am Knoten {column}')
            ax[1].legend()
            ax[1].grid()
            ax[1].set_ylabel('Strom [A]')
            ax[1].set_xlabel('Zeitpunkt [MM-TT hh]')
            ax[1].set_title('Strom über der Zeit - Ergebnisse der Optimierung')

            plt.show()



if __name__ == '__main__':
    t0 = time.time()
    test = GridLineOptimizer(6, bev_buses=list(range(6)), resolution=15)
    # print(test.buses)
    # print(test.lines)
    # print(test.impedances)
    # print(test.u_min)
    test.display_target_function()
    test.display_min_voltage_constraint()
    test.display_max_current_constraint()
    test.display_track_socs_constraint()
    res = test.run_optimization_single_timestep(tee=True)
    #test.run_optimization_rolling_horizon(24, tee=False)
    dt = time.time() - t0
    print('fertig nach', dt, 'sek')
    # test.optimization_model.I.pprint()
    # print(res)
    test.optimization_model.SOC.pprint()
    print(test.results_I)

    test.plot_results()
    test.display_ensure_final_soc_constraint()

