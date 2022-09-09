"""
Author: André Ulrich
--------------------
Klasse GridLineOptimizer, die ein Model zur Optimierung der Ladeleistungen von Ladesäulen entlang eines Netzstrahls
erzeug. Um die Ergebnisse hinterher validieren zu können, ist auch direkt ein pandapower-Netz mit enthalten, mit dem
man nach der Optimierung die Ergebnisse überprüfen kann.

Roling horizon geht jetzt!! die Ergebnisse von den SOCs schauen genauso aus, wie beim fixed horizon. Allerdings schauen
die Is komisch aus (obwohl die ja im Model mit den SOCs verknüpft sind?!) => irgendwo Fehler beim Auslesen der Is
aus dem Model??

Irgendiwe muss noch sichergestellt werden, dass die Optimierung nicht abschmiert, wenn die gewünschten SOCs zur
gewünschten Uhrzeit nicht erreicht werden könne
=> dafür erstmal die soc_lower_bounds ausschalten (bzw. den Teil, wo ab t_target dann soc_target eingetragen wird)
dann muss noch in der Zielfunktion der Anreiz geschaffen werden, dass wirklich bis zu deren t_target möglichst viel
geladen wird (momentan starten die erst kurz vor Ende) => Zielfunktion nur occupancy_times indexieren (dann muss
aber auch noch voltages entsprechend angepasst werden)

Beim upper und lower bounds von SOC und I auch mal so einstellen, dass da wirklich nur diejenigen Zeitpunkte drin sind,
an denen wirklich geladen wird **

bei prepare_soc_upper- und -lower_bounds noch dafür sorgen, dass als lower bound beim t_target vom jeweiligen BEV
auch wirklich der soc_target steht

R von 0.04 auf 0.004 (was realistischer ist, für Leitungen von ca. 15m länge und 0.255 Ohm/km spezifischem Widerstand)

Versionsgeschichte:
V.1: upper und lower bounds der Variables als dict für die einzelnen timesteps => dadurch entfällt die Subtarktion
des current_timestep beim Auslesen der bounds im model, außerdem intuitiver indexieren

V.2: "intelligentes" Set occupancy_times zum Indexieren. Darin sind jetzt nur noch diejenigen timesteps enthalten,
an denen am jeweiligen Knoten auch wirklich ein BEV steht zum Laden => dadurch kann man de if-Abfrage vor den rules
weglassen, die prüft, ob man im Ladezeitraum des entsprechenden BEVs ist. Eventuell können so auch mehrere separate
Ladevorgänge an einem Bus innerhalb eines Horizonts ermöglicht werden (weil ja dann die timesteps eigentlich nicht mehr
durchgängig miteinander verbunden sind.
"""

_pandapower_available = True
_networkx_available = True
_pandas_available = True
_matplotlib_available = True
_ipopt_available = True

import pyomo.environ as pe
import itertools as itt
from distutils.spawn import find_executable

try:
    import pandas as pd

except ModuleNotFoundError:
    print('\nWARNING: module pandas not available, some features are',
          'only available with pandas\n')
    _pandas_available = False

try:
    import pandapower as pp

except ModuleNotFoundError:
    print('\nWARNING: module pandapower not available, some features are',
          'only available with pandapower\n')
    _pandapower_available = False

try:
    import matplotlib.pyplot as plt
    import matplotlib.dates as mdates

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

if not find_executable('ipopt'):
    print('\nWARNING: ipopt solver not available, solving NLP is not',
          'supported without appropriate solver\n')
    _ipopt_available = False

import time
import os
from os import path as op
import numpy as np
from battery_electric_vehicle import BatteryElectricVehicle as BEV
from household import Household as HH



class GridLineOptimizer:
    global _pandapwer_available
    global _networkx_available
    global _pandas_available
    global _matplotlib_available
    global _ipopt_available


    _OPTIONS = {'distribute loadings': False,
                'log results': False,
                'consider linear': True,
                'fairness': 27,
                'equal SOCs': 1,
                'steady charging': (0, 0),
                'atillas constraint': False,
                'equal products': False}

    def __init__(self, number_buses, bevs, households, trafo_power, resolution, horizon_width=24,
                 voltages=None, line_impedances=None, line_lengths=None, line_capacities=None,
                 solver='glpk'):
        self.rolling = False
        self.current_timestep = 0
        self.resolution = resolution
        self.horizon_width = horizon_width
        self.number_buses = number_buses
        self.buses = self._make_buses()
        self.lines = self._make_lines()
        self._make_times()
        self.voltages = self._make_voltages(voltages)

        self.u_trafo = 400
        self.u_min = 0.91 * self.u_trafo
        self.p_trafo = trafo_power
        self.i_max = self.p_trafo * 1000 / self.u_trafo

        self.solver = solver
        self.solver_factory = pe.SolverFactory(self.solver)

        self.line_capacities = self._make_line_capacities(line_capacities)
        self.impedances = self._make_impedances(line_impedances)
        self.line_lengths = self._make_line_lengths(line_lengths)
        self.resulting_impedances = self._make_resulting_impedances()

        self._make_bev_dict(bevs)
        self._setup_bevs()
        self.households = households

        self._prepare_soc_lower_bounds()
        self._prepare_soc_upper_bounds()
        self._fix_first_socs()

        self._prepare_i_lower_bounds()
        self._prepare_i_upper_bounds()

        self._setup_model()

        # when using rolling horizon the results of the first timestep of each horizon
        # are store in here
        self.results_I = {bus: [] for bus in self.bevs}
        self.results_SOC = {bus: [] for bus in self.bevs}


    def _prepare_i_lower_bounds(self):
        """
        setup lower bounds of charging currents. They cant go below 0A.
        :return: None
        """
        i_lower_bounds = {bev.home_bus: {t: 0 for t in self.times} for bev in self.bevs.values()}
        self.i_lower_bounds = i_lower_bounds


    def _prepare_i_upper_bounds(self):
        """
        setup upper bounds of charging currents. They can only be greater than 0A
        if the according BEV is at the chargin station at the according timestep.
        They cant be greater than the charging power of the charging station
        divided by the node voltage.
        :return: None
        """
        # first, setup everything to be maximum
        i_upper_bounds = {bev.home_bus: {t: bev.p_load/0.4 for t in self.times} for bev in self.bevs.values()}
        for bev in self.bevs.values():
            # than, for each bev, set it to 0A if the bev is not at the charger at this timestep
            if self.current_timestep < bev.t_start:
                # ... all the times before t_start there has to be 0
                i_upper_bounds[bev.home_bus].update({t: 0 for t in self.times if t < bev.t_start})
                # and the same for all the times after t_target
                i_upper_bounds[bev.home_bus].update({t: 0 for t in self.times if t > bev.t_target})
            # if the current_timestep is between t_start and t_target, than ...
            elif bev.t_start <= self.current_timestep and self.current_timestep < bev.t_target:
                # ... all the times after t_target have to be 0
                i_upper_bounds[bev.home_bus].update({t: 0 for t in self.times if t > bev.t_target})
            # if the current_timestep is is above t_target ...
            elif self.current_timestep >= bev.t_start:
                # all the times have to be 0
                i_upper_bounds[bev.home_bus].update({t: 0 for t in self.times})

        self.i_upper_bounds = i_upper_bounds


    def _make_bev_dict(self, bevs):
        """
        make a nice dict out of all the BEVs, for more intuitive
        indexing (via their home_bus)
        :param bevs: list of BatteryElectricVehicles
        :return: None
        """
        bev_dict = {bev.home_bus: bev for bev in bevs}
        # just in case the BEVs have been passed in a random order
        # => sort them according to their home_bus
        self.bevs = dict(sorted(bev_dict.items()))


    def _setup_bevs(self):
        """
        tell each BEV the used width of the horizon in time
        :return: None
        """
        for bev in self.bevs.values():
            bev.set_horizon_width(self.horizon_width)
            #bev.make_occupancies()


    def _prepare_soc_lower_bounds(self):
        """
        setup lower bounds of the SOCs. They cant be less than the soc_start of the according BEV.
        :return: None
        """
        soc_lower_bounds = {bev.home_bus: {t: bev.soc_start for t in self.times} for bev in self.bevs.values()}
        self.soc_lower_bounds = soc_lower_bounds


    def _prepare_soc_upper_bounds(self):
        """
        setup upper bounds of the socs. They must ot go higher than the desired soc_target of
        the according BEV.
        :return: None
        """
        soc_upper_bounds = {bev.home_bus: {t: bev.soc_target for t in self.times} for bev in self.bevs.values()}
        self.soc_upper_bounds = soc_upper_bounds


    def _fix_first_socs(self):
        """
        make sure, that at the first considered timestep, the upper bounds equal the lower
        bounds (soc_start of the according bev). This ensures that the BEVs really start
        their charging with their soc_start (otherwise the optimizer would be allowed to
        simply raise the SOC at the beginning of charging in order to allways fullfill
        the wishes of the BEVs).
        :return: None
        """
        for bev in self.bevs.values():
            self.soc_upper_bounds[bev.home_bus][0] = bev.soc_start


    def _make_buses(self):
        """
        prepare buses list for indexing pyomo sets in the optimization model.
        :return: the desired list
        """
        return list(range(self.number_buses))


    def _make_lines(self):
        """
        prepare lines list for indexing pyomo sets in the optimization model.
        :return: the desired list
        """
        return list(range(self.number_buses))


    def _make_times(self):
        """
        prepare timesteps list for indexing pyomo sets in the optimization model.
        :return: the desired list
        """
        self.times = list(range(self.current_timestep, self.current_timestep+self.horizon_width*int(60/self.resolution)))


    def _make_line_capacities(self, capacities):
        """
        prepare line capacities dict (what current each line can conduct) for indexing
        pyomo parameters in the optimization model.
        :param capacities: list of conducting capacites (A)
        :return: the desired dict
        """
        if capacities == None:
            return {i: self.i_max for i in self.lines}

        elif type(capacities) == int or type(capacities) == float:
            return {i: capacities for i in self.lines}

        else:
            if not len(capacities) == len(self.lines):
                raise ValueError("Length of gridline is {}, but {} line capacities were passed"
                                 .format(len(self.lines), len(capacities)))

            return {i: capacities[i] for i in self.lines}


    def _make_voltages(self, voltages):
        """
        prepare node voltaged dict for indexing pyomo parameters in the optimization model
        :param voltages: list of node voltages (V)
        :return: the desired dict
        """
        if voltages == None:
            return {i: 400-(i+1)/2 for i in self.buses}

        elif type(voltages) == int or type(voltages) == float:
            return {i: voltages for i in self.buses}

        else:
            if not len(voltages) == len(self.lines):
                raise ValueError("Length of gridline is {}, but {} node voltages were passed"
                                 .format(len(self.lines), len(impedances)))

                return {i: voltages[i] for i in self.buses}



    def _make_impedances(self, impedances):
        """
        prepare specific line impedances dict for indexing pyomo parameters in the optimization model.
        :param impedances: list of specific impedances
        :return: the desired dict
        """
        if impedances == None:
            return {i: 2e-4 for i in self.lines}

        elif type(impedances) == int or type(impedances) == float:
            return {i: impedances for i in self.lines}

        else:
            if not len(impedances) == len(self.lines):
                raise ValueError("Length of gridline is {}, but {} line impedances were passed"
                                 .format(len(self.lines), len(impedances)))

            return {i: impedances[i] for i in self.lines}


    def _make_line_lengths(self, lenghts):
        """
        prepare line lengths dict for indexing pyomo parameters in the optimization model.
        :param lenghts: list of line lengths
        :return: the desired dict
        """
        if lenghts == None:
            return {i: 20 for i in self.lines}

        elif type(lenghts) == int or type(lenghts) == float:
            return {i: lenghts for i in self.lines}

        else:
            if not len(lenghts) == len(self.lines):
                raise ValueError("Length of gridline is {}, but {} line lenghts were passed"
                                 .format(len(self.lines), len(lenghts)))

            return {i: lenghts[i] for i in self.lines}


    def _make_resulting_impedances(self):
        """
        prepare resulting impedances (=line legth*line specific impedance) dict for indexing
        pyomo parameters in the optimization model.
        :return: the desired dict
        """
        return {num: self.line_lengths[num]*impedance for num, impedance in enumerate(self.impedances.values())}


    def _carry_over_last_socs(self):
        """
        after each run of the optimizer, fetch the results from the optimization
        model (the socs of the second timestep in the current horizon) and use
        these values as new upper and lower bounds for the first timestep in
        the next horizon ("sandwich method" to ensure energy conservation over
        multiple horizons). IMPORTANT: call it AFTER the current_timestep has
        been incremented and AFTER _prepare_soc_upper/lower_bounds have been
        called
        :return: None
        """
        # get the socs out of the optimization model
        socs_to_carry_over = [self.optimization_model.SOC[self.current_timestep,
                              bus].value for bus in self.bevs]

        # and use them as upper and lower bounds, first lower
        for num, bus in enumerate(self.bevs):
            self.soc_lower_bounds[bus][self.current_timestep] = socs_to_carry_over[num]
            self.soc_upper_bounds[bus][self.current_timestep] = socs_to_carry_over[num]


        self._update_bev_socs(socs_to_carry_over)


    def _prepare_next_timestep(self):
        """
        when using the rolling horizon, make sure some things are well prepared
        befor start building the optimization model for the next horizon:
        1. increment the current_timestep
        2. build according times list for correct indexing in the next horizons
        optimization model
        3. prepare all the variables upper and lower bounds accordingly
        4. for ensuring energy conservation over multiple horizos: take
        the socs of the second timestep of the current horizon and use
        them as upper and lower bounds ("sandwich method") for the
        first timestep soc in the next horizon.
        :return: None
        """
        if not self.rolling:
            self.rolling = True

        self.current_timestep += 1
        self._make_times()

        self._prepare_soc_lower_bounds()
        self._prepare_soc_upper_bounds()

        self._carry_over_last_socs()

        self._prepare_i_lower_bounds()
        self._prepare_i_upper_bounds()


    def _update_bev_socs(self, values):
        """
        Assigns the SOCs to the corresponding BEVs instances
        :param values: List with SOCs for all the BEVs
        :return: None
        """
        for num, bev in enumerate(self.bevs.values()):
            bev.update_soc(values[num])


    def _setup_model(self):
        """
        create optimization model for the considered time horizon.
        1. build sets for indexing all the other components
        2. build parameters, indexed in the according sets
        3. build decission variables, indexed in the according sets
        4. build restrictions, indexed in the according sets
        :return: prpepared optimization model instance
        """
        model = pe.ConcreteModel('GridLineOptimization')

        # create
        model.buses = pe.Set(initialize=self.buses)
        model.lines = pe.Set(initialize=self.lines)
        model.times = pe.Set(initialize=self.times)
        model.charger_buses = pe.Set(initialize=[bev.home_bus for bev in self.bevs.values()])

        # create parameters
        model.impedances = pe.Param(model.lines, initialize=self.resulting_impedances)
        model.voltages = pe.Param(model.buses, initialize=self.voltages)
        model.u_min = self.u_min
        model.u_trafo = self.u_trafo
        model.i_max = self.i_max
        model.line_capacities = pe.Param(model.lines, initialize=self.line_capacities)

        def get_household_currents(model, time, bus):
            """
            get  household currents, gets passed to
            constructor of pyomo Parameter. Must not be called on its own.
            :param model:
            :param time: 
            :param bus: 
            :return: household currents
            """""
            return self.households[bus].load_profile[time] / self.voltages[bus]


        model.household_currents = pe.Param(model.times*model.buses, initialize=get_household_currents,
                                            mutable=True)

        def get_soc_bounds(model, time, bus):
            """
            get upper/lower bounds for charging currents, gets passed to
            constructor of pyomo Variable. Must not be called on its own.
            :param model:
            :param time:
            :param bus:
            :return: upper/lower bounds
            """
            return (self.soc_lower_bounds[bus][time], self.soc_upper_bounds[bus][time])


        def get_i_bounds(model, time, bus):
            """
            get upper/lower bounds for SOCs, gets passed to
            constructor of pyomo Variable. Must not be called on its own.
            :param model:
            :param time:
            :param bus:
            :return: upper/lower bounds
            """
            return (self.i_lower_bounds[bus][time], self.i_upper_bounds[bus][time])


        # create decission variables
        model.I = pe.Var(model.times*model.charger_buses, domain=pe.NonNegativeReals, bounds=get_i_bounds)
        model.SOC = pe.Var(model.times*model.charger_buses, domain=pe.PositiveReals, bounds=get_soc_bounds)

        # create target function
        def max_power_rule(model):
            """
            create expression for target function. Gets passed to constructor of pyomo Objective.
            Must not be called on its own.
            :param model:
            :return: the expression
            """
            return sum(sum(model.I[t, b] for t in model.times) for b in model.charger_buses)


        model.max_power = pe.Objective(rule=max_power_rule, sense=pe.maximize) # maximize the currents

        # create restrictions
        def min_voltage_rule(model, t):
            """
            create expression for constraint: the voltage at the last node must not fall
            below tolerable voltage band. Gets passed to conrtructor of pyomo Constraint.
            Must not be caled on its own.
            :param model:
            :param t:
            :return: the expression
            """
            # lambda functions just for accessing the l counter inside sum, this way we can
            # avoid one very long expression.
            hh_currs = lambda l: sum(model.household_currents[t, n] for n in model.buses if n >= l)
            bev_currs = lambda l: sum(model.I[t, n] for n in model.charger_buses if n >= l)
            # the actual expression, more compact
            return model.u_trafo - sum(model.impedances[l] * (hh_currs(l) + bev_currs(l))
                                       for l in model.lines) >= model.u_min


        def max_current_rule(model, t):
            """
            create expression for constraint: the sum of all currents must not exceed the
            maximum current the transformer can conduct (P_trafo/u_0). Gets passed to conrtructor
            of pyomo Constraint. Must not be caled on its own.
            :param model:
            :param t:
            :return: the expression
            """
            return sum(model.I[t, b] for b in model.charger_buses) +\
                   sum(model.household_currents[t, b] for b in model.buses) <= model.i_max


        def line_capacities_rule(model, t, b):
            fn = b # first relevant node to consider
            return sum(model.I[t, b] for b in model.charger_buses if b >= fn) + sum(
                model.household_currents[t, b] for b in model.buses if b >= fn
            ) <= model.line_capacities[b]


        def track_socs_rule(model, t, b):
            # schauen, dass man immer nur bis zum vorletzten timestep geht (weil es
            # sonst kein t+1 mehr geben würde beim letzten timestep)
            if t < self.current_timestep + self.horizon_width * int(60/self.resolution) - 1:#23:
            #if t >= self.bevs[b].t_start and t <= self.bevs[b].t_target:
                return (model.SOC[t, b] + model.I[t, b] * model.voltages[b] * self.resolution/60 / 1000
                        / self.bevs[b].e_bat*100 - model.SOC[t+1, b]) == 0

            else:
                return pe.Constraint.Skip


        def fair_charging_rule(model, t, b, pb):
            if b > pb:
                return model.I[t, b] - model.I[t, pb] <=\
                       type(self)._OPTIONS['fairness']

            else:
                return pe.Constraint.Skip


        def equal_socs_rule(model, j, k):
            if self.rolling:
                ft_j = self.bevs[j].t_target#self.current_timestep + self.horizon_width * 60 / self.resolution - 1
                ft_k = self.bevs[k].t_target
                #if j < np.max(model.charger_buses) + 1: #hier muss nicht die Länge, sondern das MAXIMUM in charger_buses verwendet werden!
                if self.current_timestep < min(ft_j, ft_k):
                    if j > k:
                        fullfillment_j = (model.SOC[ft_j, j] - self.bevs[j].soc_start)/(self.bevs[j].soc_target-self.bevs[j].soc_start)
                        fullfillment_k = (model.SOC[ft_k, k] - self.bevs[k].soc_start)/(self.bevs[k].soc_target-self.bevs[k].soc_start)
                        return fullfillment_j - fullfillment_k <= type(self)._OPTIONS['equal SOCs']
                    else:
                        return pe.Constraint.Skip

                else:
                    return pe.Constraint.Skip

            else:
                ft = self.current_timestep + self.horizon_width * 60 / self.resolution - 1
                if j > k:
                    fullfillment_j = (model.SOC[ft, j] - self.bevs[j].soc_start) / (self.bevs[j].soc_target - self.bevs[j].soc_start)
                    fullfillment_k = (model.SOC[ft, k] - self.bevs[k].soc_start) / (self.bevs[k].soc_target - self.bevs[k].soc_start)
                    return fullfillment_j - fullfillment_k <= type(self)._OPTIONS['equal SOCs']
                else:
                    return pe.Constraint.Skip


        def equal_product_rule(model, b):
            """
            the product of bev soc_start and the amount of what they reach
            of their intended loaded difference, should be equal for all bevs
            => people with lower start_socs get more loaded (relative to how much
            they want to be loaded
            :param model:
            :param b:
            :return:
            """
            lb = model.charger_buses.prevw(b)
            # last considered timestep
            lt_b = self.bevs[b].t_target#self.current_timestep + self.horizon_width * 60/self.resolution -1
            lt_lb = self.bevs[lb].t_target
            # first considered timestep
            #ft = self.bevs[b].t_start
            if  self.current_timestep < min(self.bevs[b].t_target, self.bevs[lb].t_target):
                if b > 0:
                    fullfillment_b = (model.SOC[lt_b, b] - self.bevs[b].soc_start) / (self.bevs[b].soc_target - self.bevs[b].soc_start)
                    fullfillment_lb = (model.SOC[lt_lb, lb] - self.bevs[lb].soc_start) / (self.bevs[lb].soc_target - self.bevs[lb].soc_start)
                    return self.bevs[b].soc_start * fullfillment_b - self.bevs[lb].soc_start * fullfillment_lb == 0

                else:
                    return pe.Constraint.Skip

            else:
                return pe.Constraint.Skip


        def atillas_rule(model, b):
            """
            function to be passed to pyomo Constructor for Constraints - not intended
            to be called on its own.
            :param model: model reference
            :param b: cycles through model.charger_buses
            :param t: cycles through model.times
            :return: expression for constraint
            """
            lb = model.charger_buses.prevw(b)
            ct = self.current_timestep # current timestep
            #if t >= self.bevs[b].t_start and t < self.bevs[b].t_target:
            if ct < self.bevs[b].t_target:
                if b > 0:
                    ts_lb = self.bevs[lb].t_start
                    tt_lb = self.bevs[lb].t_target
                    currents_lb = sum(model.I[t, lb] for t in model.times if t >= ts_lb and t < tt_lb)
                    dt_lb = self.bevs[lb].t_target - ct
                    dsoc_lb = self.bevs[lb].soc_target - self.bevs[lb].current_soc

                    ts_b = self.bevs[b].t_start
                    tt_b = self.bevs[b].t_target
                    currents_b = sum(model.I[t, b] for t in model.times if t >= ts_b and t < tt_b)
                    dt_b = self.bevs[b].t_target - ct
                    dsoc_b = self.bevs[b].soc_target - self.bevs[b].current_soc
                    # return sum(model.I[t, lb] for t in model.times if t >= self.bevs[lb].t_start and t < self.bevs[lb].t_target) * \
                    #        ((self.bevs[lb].t_target - ct) / (self.bevs[lb].soc_target - self.bevs[lb].current_soc))**2 \
                    # <= sum(model.I[t, b] for t in model.times if t >= self.bevs[lb].t_start and t < self.bevs[lb].t_target) * ((self.bevs[b].t_target - ct) / \
                    #                                              (self.bevs[b].soc_target - self.bevs[b].current_soc))**2
                    return currents_lb * (dt_lb / dsoc_lb)**2 <= currents_b * (dt_b / dsoc_b)**2

                else:
                    return pe.Constraint.Skip

            else:
                return pe.Constraint.Skip


        model.min_voltage = pe.Constraint(model.times, rule=min_voltage_rule)
        model.max_current = pe.Constraint(model.times, rule=max_current_rule)
        model.keep_line_capacities = pe.Constraint(model.times*model.lines, rule=line_capacities_rule)
        model.track_socs = pe.Constraint(model.times*model.charger_buses, rule=track_socs_rule)

        if type(self)._OPTIONS['fairness'] < 27:
            model.fair_charging = pe.Constraint(model.times*model.charger_buses*model.charger_buses, rule=fair_charging_rule)
            #model.fair_charging.pprint()

        if type(self)._OPTIONS['equal SOCs'] < 1:
            model.equal_socs = pe.Constraint(model.charger_buses*model.charger_buses, rule=equal_socs_rule)
            #model.equal_socs.pprint()

        if type(self)._OPTIONS['atillas constraint'] == True:
            model.atillas_constraint = pe.Constraint(model.charger_buses, rule=atillas_rule)
            #model.atillas_constraint.pprint()

        if type(self)._OPTIONS['equal products'] == True:
            model.equal_products = pe.Constraint(model.charger_buses, rule=equal_product_rule)

        #return model
        self.optimization_model = model


    def _setup_grid(self):
        """
        create pandapower grid to check the results of optimization. Not needed anymore, because
        the EMO is used for simulating the grid
        :return: pandapower grid
        """
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
            pp.create_transformer_from_parameters(grid, hv_bus=0, lv_bus=1, sn_mva=self.p_trafo/1000,
                                                  vn_hv_kv=20, vn_lv_kv=0.4, vkr_percent=1.5, pfe_kw=0.4,
                                                  i0_percent=0.4, vk_percent=6)

            # Leitungen erzeugen
            for nr in range(self.number_buses-1):
                pp.create_line_from_parameters(grid, from_bus=nr+1, to_bus=nr+2, r_ohm_per_km=1,
                                               length_km=1/self.impedances[nr], name='line '+str(nr),
                                               x_ohm_per_km=0, c_nf_per_km=0, max_i_ka=0.142)

            return grid


    def export_grid(self, filename):
        """
        create an excel file called 'optimized_grid.xlsx' containing all the
        information of the optimized grid, that the EMO needs to construct a
        pandas grid from, using LowVoltageSystem.make_system_from_excel_file
        :return: None
        """
        num_buses = 2 + self.number_buses
        # for the sheet 'Lines'
        line_no = [i for i in range(num_buses-2)]
        from_bus = [i+1 for i in range(num_buses-2)]
        to_bus = [i+2 for i in range(num_buses-2)]
        length = [15 for _ in range(num_buses-2)]

        lines_dict = {'Line No.': line_no,
                      'From Bus': from_bus,
                      'To Bus': to_bus,
                      'Length': length}

        lines_df = pd.DataFrame(lines_dict)

        # for the sheet 'Buses'
        home_buses = [i+2 for i in self.bevs.keys()] # helper, since enumeration for bevs begins at 0
        bus_no = [i for i in range(num_buses)]
        x = bus_no # +1 in x-direction for each bus
        y = [0 for _ in range(num_buses)]
        household = ['No' if i < 2 else 'Yes' for i in range(num_buses)]
        wallbox = ['Yes' if i in home_buses else 'No' for i in range(num_buses)]

        buses_dict = {'Bus No.': bus_no,
                      'X': x,
                      'Y': y,
                      'Household': household,
                      'Wallbox': wallbox}

        buses_df = pd.DataFrame(buses_dict)

        # write to excel file
        with pd.ExcelWriter('{}.xlsx'.format(filename)) as writer:
            lines_df.to_excel(writer, sheet_name='Lines', index=False)
            buses_df.to_excel(writer, sheet_name='Busses', index=False)


    def export_household_profiles(self):
        num_timesteps = int(self.horizon_width * 60 / self.resolution)
        return {household.home_bus: household.load_profile[0:num_timesteps] for household in self.households}


    def export_I_results(self):
        if not self.rolling:
            return {bev: [self.optimization_model.I[t, bev].value for t in self.times]
                    for bev in self.bevs.keys()}

        else:
            return self.results_I


    def get_grid_specs(self):
        specs = {'buses': self.number_buses,
                 'S transformer': self.p_trafo,
                 'line specific impedances': self.impedances,
                 'line lenghts': self.line_lengths,
                 'line resulting impedances': self.resulting_impedances,
                 'line capacities': self.line_capacities}

        return specs


    def display_target_function(self):
        self.optimization_model.max_power.pprint()


    def display_min_voltage_constraint(self):
        self.optimization_model.min_voltage.pprint()


    def display_max_current_constraint(self):
        self.optimization_model.max_current.pprint()


    def display_keep_line_capacities_constraint(self):
        self.optimization_model.keep_line_capacities.pprint()


    def display_track_socs_constraint(self):
        self.optimization_model.track_socs.pprint()


    def display_ensure_final_soc_constraint(self):
        self.optimization_model.ensure_final_soc.pprint()


    def log_results(self):
        """
        logs the results of the optimization to an external csv-file.
        :return: None
        """
        # first check, if csv already exists
        if not op.isdir('../Data/Results'):
            os.mkdir('../Data/Results')

        # create dict with optimization results
        # results = {'timestep': self.current_timestep,
        #            'SOC [%]': [self.optimization_model.SOC[self.current_timestep, bev].value for bev in self.bevs],
        #            'I [A]': [self.optimization_model.I[self.current_timestep, bev].value for bev in self.bevs]
        #            }

        results = {f'current at node {bev}': [self.optimization_model.I[self.current_timestep+1, bev].value]
                   for bev in self.bevs}
        results['timestep'] = [self.current_timestep]

        results = pd.DataFrame(results)

        if not op.isfile('../Data/Results/results.csv'):
            results.to_csv('../Data/Results/results.csv', index=False)

        else:
            results.to_csv('../Data/Results/results.csv', index=False, header=False, mode='a')



    # nach jedem Optimierungsdurchlauf die Ergebnisse aus dem Model und
    # die SOCs in den BEVs speichern, I hier irgendwo...
    def _store_results(self):
        """
        Fragt die Optimierungsergebnisse der Entscheidungsvariablen I und SOC aus dem
        Optimierungsmodel ab und speichert diese.
        :return:
        """
        for bus in self.bevs:  # das liefert ja die home_buses
            # Werte aus model abfragen
            SOC = self.optimization_model.SOC[self.current_timestep, bus].value # +1 dazu
            I = self.optimization_model.I[self.current_timestep, bus].value # +1 dazu
            # und in Ergebnisliste eintragen
            self.results_SOC[bus].append(SOC)
            self.results_I[bus].append(I)


    # simuliert einen Tag mit den Werten für einen Tag, fixer Horizont
    def run_optimization_single_timestep(self, **kwargs):
        self.solver_factory.solve(self.optimization_model, tee=kwargs['tee'])
        if type(self)._OPTIONS['log results']:
            self.log_results()


    def run_optimization_rolling_horizon(self, complete_horizon, **kwargs):
        steps = int(complete_horizon * 60/self.resolution)
        for i in range(steps):
            print(i)
            self.run_optimization_single_timestep(tee=kwargs['tee'])
            self._store_results()
            self._prepare_next_timestep()#=kwargs['update_bevs'])
            self._setup_model()


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


    def plot_I_results(self, legend=True, save=False, usetex=False,
                       compact_x=False, **kwargs):
        if usetex:
            plt.rcParams['text.usetex'] = True
            plt.rcParams['font.family'] = 'serif'
            plt.rcParams['grid.linewidth'] = 0.4
            plt.rcParams['lines.linewidth'] = 1
            plt.rcParams['legend.fontsize'] = 8
            plt.rcParams['font.size'] = 11

        if compact_x:
            x_fmt = mdates.DateFormatter('%H')

        Is = {bus: [] for bus in self.bevs}
        for time in self.times:
            for bus in self.bevs:
                Is[bus].append(self.optimization_model.I[time, bus].value)

        Is_df = pd.DataFrame(Is)
        Is_df.index = pd.date_range(start='2021', periods=len(Is_df), freq=str(self.resolution) + 'min')

        fig, ax = plt.subplots(1, 1, figsize=(6.5, 1.75))

        for column in Is_df.columns:
            ax.plot(Is_df.index, Is_df[column], marker=kwargs['marker'], label=f'Strom zum BEV am Knoten {column}')
        if legend:
            ax.legend()
        ax.grid()
        ax.set_ylabel('Ladestrom [A]')
        if compact_x:
            ax.xaxis.set_major_formatter(x_fmt)
            ax.set_xlabel('Zeit [hh]')
        else:
            ax.set_xlabel('Zeitpunkt [MM-TT hh]')
        #ax.set_title('Strom über der Zeit -- Ergebnis der Optimierung')
        if not save:
            plt.show()

        else:
            plt.savefig('res_opt_i.pdf', bbox_inches='tight')


    def plot_SOC_results(self, legend=True, save=False, usetex=False,
                         compact_x=False, **kwargs):
        if usetex:
            plt.rcParams['text.usetex'] = True
            plt.rcParams['font.family'] = 'serif'
            plt.rcParams['grid.linewidth'] = 0.4
            plt.rcParams['lines.linewidth'] = 1
            plt.rcParams['legend.fontsize'] = 8
            plt.rcParams['font.size'] = 11

        if compact_x:
            x_fmt = mdates.DateFormatter('%H')

        SOCs = {bus: [] for bus in self.bevs}
        for time in self.times:
            for bus in self.bevs:
                SOCs[bus].append(self.optimization_model.SOC[time, bus].value)

        SOCs_df = pd.DataFrame(SOCs)
        SOCs_df.index = pd.date_range(start='2021', periods=len(SOCs_df), freq=str(self.resolution) + 'min')

        fig, ax = plt.subplots(1, 1, figsize=(6.5, 1.75))

        for column in SOCs_df.columns:
            ax.plot(SOCs_df.index, SOCs_df[column], marker=kwargs['marker'], label=f'SOC des BEV am Knoten {column}')
        if legend:
            ax.legend()
        ax.grid()
        if usetex:
            ax.set_ylabel('SOC [\%]')
        else:
            ax.set_ylabel('SOC [%]')

        if compact_x:
            ax.xaxis.set_major_formatter(x_fmt)
            ax.set_xlabel('Time [hh]')
        else:
            ax.set_xlabel('Zeitpunkt [MM-TT hh]')
        #ax.set_title('SOC über der Zeit -- Ergebnis der Optimierung')
        if not save:
            plt.show()

        else:
            plt.savefig('res_opt_soc.pdf', bbox_inches='tight')


    def plot_all_results(self, legend=True, save=False, usetex=False, compact_x=False,
                         export_data=False, **kwargs):
        if usetex:
            plt.rcParams['text.usetex'] = True
            plt.rcParams['font.family'] = 'serif'
            plt.rcParams['grid.linewidth'] = 0.4
            plt.rcParams['lines.linewidth'] = 1
            plt.rcParams['legend.fontsize'] = 8
            plt.rcParams['font.size'] = 10.95

        # if compact_x:
        #     x_fmt = mdates.DateFormatter('%H')

        if not _pandas_available or not _matplotlib_available:
            print('\nWARNING: unable to plot results\n')

        else:
            Is_df, SOCs_df = self._gather_data_for_plotting()
            fig, ax = plt.subplots(2, 1, figsize=(6.3, 4), sharex=True)
            for column in SOCs_df.columns:
                ax[0].plot(SOCs_df.index, SOCs_df[column], marker=kwargs['marker'], label=f'Knoten {column+1}')
            if legend:
                ax[0].legend()
            ax[0].grid()
            if usetex:
                ax[0].set_ylabel('SOC [\%]')
            else:
                ax[0].set_ylabel('SOC [%]')
            #ax[0].set_title('SOC over time - results of optimization', fontsize=20)

            for column in Is_df.columns:
                ax[1].plot(Is_df.index, Is_df[column], marker=kwargs['marker'], label=f'Knoten {column+1}')
            if legend:
                ax[1].legend()
            ax[1].grid()
            ax[1].set_ylabel('Strom [A]')
            ax[1].set_xlabel('Time [mm-dd hh]', fontsize=11)
            if compact_x:
                x_fmt = mdates.DateFormatter('%H')
                ax[1].xaxis.set_major_formatter(x_fmt)
                ax[1].set_xlabel('Zeit [hh]')
            #ax[1].set_title('Current over time - results of optimization', fontsize=20)

            if export_data:
                SOCs_df.index = np.linspace(0, self.horizon_width, self.horizon_width * int(60 / self.resolution))
                SOCs_df.to_csv('SOCs.dat', sep='\t')
                Is_df.index = np.linspace(0, self.horizon_width, self.horizon_width * int(60 / self.resolution))
                Is_df.to_csv('Is.dat', sep='\t')

            if not save:
                plt.show()

            else:
                plt.savefig('res_opt.pdf', bbox_inches='tight')


    def _gather_data_for_plotting(self):
        if not self.rolling:
            # if there was no rolling horizon, get results out of optimization model
            SOCs = {bus: [self.optimization_model.SOC[time, bus].value for time in self.times] for bus in self.bevs}
            Is = {bus: [self.optimization_model.I[time, bus].value for time in self.times] for bus in self.bevs}

            SOCs_df = pd.DataFrame(SOCs)
            Is_df = pd.DataFrame(Is)
            SOCs_df.index = pd.date_range(start='2021', periods=len(SOCs_df), freq=str(self.resolution) + 'min')
            Is_df.index = pd.date_range(start='2021', periods=len(Is_df), freq=str(self.resolution) + 'min')

        else:
            # if there was the rolling horizon, get the data out of self.results_...
            SOCs_df = pd.DataFrame(self.results_SOC)
            Is_df = pd.DataFrame(self.results_I)
            SOCs_df.index = pd.date_range(start='2021', periods=len(SOCs_df), freq=str(self.resolution) + 'min')
            Is_df.index = pd.date_range(start='2021', periods=len(Is_df), freq=str(self.resolution) + 'min')

        return Is_df, SOCs_df


    def _get_socs_fullfillment(self, optimized):
        final_timestep = self.horizon_width * int(60 / self.resolution) - 1
        if optimized:
            if not self.rolling:
                final_socs = {bev: [(self.optimization_model.SOC[final_timestep, bev].value - self.bevs[bev].soc_start) / (self.bevs[bev].soc_target - self.bevs[bev].soc_start) * 100] for bev in self.bevs}
                return final_socs

            else:
                final_socs = {bev.home_bus: [(bev.current_soc - bev.soc_start) / (bev.soc_target - bev.soc_start) * 100] for bev in self.bevs.values()}
                return final_socs

        else:
            final_socs = {bev.home_bus: [(bev.current_soc - bev.soc_start) / (bev.soc_target - bev.soc_start) * 100] for bev in self.bevs.values()}
            return final_socs


    def export_socs_fullfillment(self, optimized=True):
        final_socs = self._get_socs_fullfillment(optimized)
        final_socs = pd.DataFrame(final_socs).T
        final_socs.index += 1
        final_socs.to_csv('socs_fullfillment.dat', sep='\t', index=True, header=False)


    def resample_data(func):
        def inner(self, *args, **kwargs):
            data = func(self, *args, **kwargs)
            data = pd.DataFrame(data)
            # datetimeIndex for resampling
            data.index = pd.date_range(start='2020', periods=len(data), freq=str(self.resolution)+'min')
            # do the resampling
            data_res = data.resample('1min').interpolate()
            return data_res

        return inner


    @resample_data
    # deprecated, use export_I_results instead
    def provide_data(self, dtype):
        """
        provides data (the currents for each charger at each timestep) in a fashion that is digestable for
        the EMO.sim_handler.run_sim => dict of dicts or dict of lists
        e.g.:
        {                                        or: {
        bus1: {ts1: x, ts2: x, ...., tsn: x},         bus1: [x, x, ..., x],
        bus2: {ts1: x, ts2: x, ..., tsn: x},          bus2: [x, x, ..., x],
         .                                             .
         .                                             .
         .                                             .
        busn: {ts1: x, ts2: x, ..., tsn: x}           busn: [x, x, ..., x]
        }                                            }
        :return: dict
        """
        if dtype not in ['dict', 'list']:
            raise ValueError("dtype needs to be 'dict' or 'list'.")

        if dtype == 'dict':
            return {bus: {time: self.optimization_model.I[time, bus].value for time in self.times}
                    for bus in self.optimization_model.charger_buses}

        elif dtype == 'list':
            return {bus: [self.optimization_model.I[time, bus].value for time in self.times]
                    for bus in self.optimization_model.charger_buses}


    def export_current_I_results(self, width):
        ts = self.current_timestep
        return {node: [self.optimization_model.I[time, node].value for time in range(ts, ts+width)] for node in self.bevs}


    @classmethod
    def set_options(cls, key, value):
        cls._OPTIONS[key] = value




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

