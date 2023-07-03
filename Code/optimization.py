"""
Author: André Ulrich
--------------------
Class GridLineOptimizer: All the functionality to build an optimization model according
to a specified grid line (number of buses, trafo power, number of BEV charging points,
line impedances etc.) and solve for timeseries of optimized BEV charging currents and
SOCs.

Usage: first create BatteryElectricVehicle and Household instances as needed and further
specify the grid topology (number of nodes, line impedances, trafo power...).
The correspondig optimization model gets automatically created on creation of the
GridLineOptimizer instance. Before creation, some options may be passed via class method,
to determine how the optimization model is going to look (further additional restrictions
for fair charging).
Once the instance is created, call run_optimization_fixed_horizon (for fixed horizon) or
run_optimization_rolling_horizon (for rolling horizon) to solve the optimization model.
With plot_all_results the results of optimization can be plotted and also exported.
With export_soc_fullfillments the SOC fullfillments of all the BEVs can be exported.
With export_I_results the results of optimized BEV charging currents can be exported
for further usage with EMO grid simulation.
With export_household_profiles the household load profiles can be exported for further
usage with EMO grid simulation.
With export_grid and export_grid_specs the whole grid topology can be exported for
further usage with EMO grid simulation.


Version history (only the most relevant points, full history is available on github):
-------------------------------------------------------------------------------------------------
V.1: first working formulation of optimization model

V.2: further methods for exporting grid topology etc. for the interface to EMO grid simulation

V.3: added further constraints for fair charging

V.4: this whole documentation added and doctsrings for all the (relevant) methods and some more
explanatory comments.

all the other commits in much more detail are available here:
https://github.com/AndreGismo/Masterarbeit/tree/submission)
"""

# these variables are used to tell the class wheather the according modules are available.
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
    import matplotlib.animation
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


    _OPTIONS = {'log results': False,
                'fairness': 27,
                'equal SOCs': 1,
                'atillas constraint': False,
                'equal products': False}

    def __init__(self, number_buses, bevs, households, trafo_power, resolution, horizon_width=24,
                 voltages=None, line_impedances=None, line_lengths=None, line_capacities=None,
                 solver='glpk'):
        """
        create GridLineOptimizer

        :param number_buses: number of buses in the grid line
        :param bevs: list of BatteryElectricVehivle instances
        :param households: list of Household instances
        :param trafo_power: power of transformer [kW]
        :param resolution: resolution in time [minutes]
        :param horizon_width: width of the optimized horizon [hours]
        :param voltages: voltages at the buses [V]
        :param line_impedances: specific impedances of the lines [ohm*m-1]
        :param line_lengths: length of lines [m]
        :param line_capacities: maximum current the lines can conduct [A]
        :param solver: solver to use (currently only 'glpk' makes sense)
        """
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
        self.u_min = 0.91 * self.u_trafo #0.945
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

        :param capacities: list of conducting capacities (A)
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
        prepare node voltaged dict for indexing pyomo parameters in the optimization model.

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
        called.

        :return: None
        """
        # get the socs out of the optimization model (its directly at the second
        # timestep, because current_timestep has been already incremented before
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
        Assigns the SOCs to the corresponding BEVs instances.

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

        :return: None
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
            maximum current the transformer can conduct (P_trafo/u_0). Gets passed to constructor
            of pyomo Constraint. Must not be caled on its own.

            :param model:
            :param t:
            :return: the expression
            """
            return sum(model.I[t, b] for b in model.charger_buses) +\
                   sum(model.household_currents[t, b] for b in model.buses) <= model.i_max


        def line_capacities_rule(model, t, b):
            """
            create expression for constraint: the sum of all currents flowing through
            one specific line mus not exceed this lines current capacity. Gets passed
            to constructor of pyomo Constraint. Must not be called on its own.

            :param model:
            :param t:
            :param b:
            :return: the expression
            """
            fn = b # first relevant node to consider
            return sum(model.I[t, b] for b in model.charger_buses if b >= fn) + sum(
                model.household_currents[t, b] for b in model.buses if b >= fn
            ) <= model.line_capacities[b]


        def track_socs_rule(model, t, b):
            """
            create expression for constraint: ensure energy conservation while
            charging. Gets passed to constructor of pyomo Constraint. Must not
            be called on its own.

            :param model:
            :param t:
            :param b:  the expression
            :return:
            """
            # only go till one timestep befor the final one, because otherwise, there
            # would be no next timestep.
            if t < self.current_timestep + self.horizon_width * int(60/self.resolution) - 1:
                return (model.SOC[t, b] + model.I[t, b] * model.voltages[b] * self.resolution / 60 / 1000
                        / self.bevs[b].e_bat * 100 - model.SOC[t+1, b]) == 0

            else:
                return pe.Constraint.Skip


        def fair_charging_rule(model, t, b, pb):
            """
            create expression for constraint: ensure that for all considered timesteps
            the maximum difference between the BEVs charging currents is smaller than
            a permittable maximum difference (determined by the option 'fairness').
            Gets passed to constructor of pyomo Constraint. Must not be called on
            its own.

            :param model:
            :param t:
            :param b:
            :param pb:
            :return: the expression
            """
            # smart pairwise comparison: no need to compare BEV 2 with BEV 1, if before
            # there has already been the comparison of BEV 1 with BEV 2 etc...
            if b > pb:
                return model.I[t, b] - model.I[t, pb] <=\
                       type(self)._OPTIONS['fairness']

            else:
                return pe.Constraint.Skip


        def equal_socs_rule(model, j, k):
            """
            create expression for constraint: the target fullfillment (the amount of energy
            a BEV is actually charged divided by the amount of energy it wishes to be charged)
            between the BEVs must not deviate more than a maximum permittable difference
            (determined by option 'equal socs'). Gets passed to constructor of pyomo Constraint.
            Must not be called on its own.

            :param model:
            :param j:
            :param k:
            :return: the expression
            """
            if self.rolling:
                ft_j = self.bevs[j].t_target
                ft_k = self.bevs[k].t_target
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
            create expression for constraint: the lower the start_soc of a BEV, the higher
            target fullfillment this BEV is allowed to get. Gets passed to constructor of
            pyomo Constraint. Must not be called on its own.

            :param model:
            :param b:
            :return: the expression
            """
            lb = model.charger_buses.prevw(b)
            lt_b = self.bevs[b].t_target
            lt_lb = self.bevs[lb].t_target
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
            create expression for constraint: the longer a BEV is already waiting to be charged,
            the more it should be carged. This comes from Atilla Koese bachelor thesis and is
            adapted to work with pyomo. Gets passed to constructor of pyomo Constraint. Must
            not be called on its own

            :param model:
            :param b:
            :return: the expression
            """
            lb = model.charger_buses.prevw(b) # the prevw bus
            ct = self.current_timestep # current timestep
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

                    return currents_lb * (dt_lb / dsoc_lb)**2 <= currents_b * (dt_b / dsoc_b)**2

                else:
                    return pe.Constraint.Skip

            else:
                return pe.Constraint.Skip

        # create the actual constraints
        model.min_voltage = pe.Constraint(model.times, rule=min_voltage_rule)
        model.max_current = pe.Constraint(model.times, rule=max_current_rule)
        model.keep_line_capacities = pe.Constraint(model.times*model.lines, rule=line_capacities_rule)
        model.track_socs = pe.Constraint(model.times*model.charger_buses, rule=track_socs_rule)

        # and here the optional constraints
        if type(self)._OPTIONS['fairness'] < 27:
            model.fair_charging = pe.Constraint(model.times*model.charger_buses*model.charger_buses, rule=fair_charging_rule)

        if type(self)._OPTIONS['equal SOCs'] < 1:
            model.equal_socs = pe.Constraint(model.charger_buses*model.charger_buses, rule=equal_socs_rule)

        if type(self)._OPTIONS['atillas constraint'] == True:
            model.atillas_constraint = pe.Constraint(model.charger_buses, rule=atillas_rule)

        if type(self)._OPTIONS['equal products'] == True:
            model.equal_products = pe.Constraint(model.charger_buses, rule=equal_product_rule)

        self.optimization_model = model


    def _setup_grid(self):
        """
        create pandapower grid to check the results of optimization. Not needed anymore, because
        the EMO is used for simulating the grid.

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
        pandas grid from, using LowVoltageSystem.grid_from_GLO.

        :return: None
        """
        # two more buses for trafo lv and mv/slack
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
        """
        export the load profile of the households for further use in grid simulation
        with EMO.

        :return: dict of household power demand
        """
        num_timesteps = int(self.horizon_width * 60 / self.resolution)
        return {household.home_bus: household.load_profile[0:num_timesteps] for household in self.households}


    def export_I_results(self):
        """
        export the results of the optimized BEV charging currents for further usage
        in EMO grid simulation.

        :return: dict of optimized BEV charging currents
        """
        if not self.rolling:
            # if using a fixed horizon, just get results from the optimization model
            # Variable instances
            return {bev: [self.optimization_model.I[t, bev].value for t in self.times]
                    for bev in self.bevs.keys()}

        else:
            # if using rolling horizon get the results from results_I
            return self.results_I


    def get_SOC_results(self):
        """
        get the results of SOCs of all EVs for the current horizon
        :param self:
        :return: dict of EV SOCs
        """
        return {
            bev: [
                self.optimization_model.SOC[t, bev].value
                for t in self.times
            ]
            for bev in self.bevs
        }


    def get_grid_specs(self):
        """
        export further specs of the grid for further usage with
        Low_Voltage_System.grid_from_GLO.

        :return: dict
        """
        specs = {'buses': self.number_buses,
                 'S transformer': self.p_trafo,
                 'line specific impedances': self.impedances,
                 'line lenghts': self.line_lengths,
                 'line resulting impedances': self.resulting_impedances,
                 'line capacities': self.line_capacities}

        return specs


#### these can be used for debugging purposes: do the constraints look as desired? ####################
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

#######################################################################################################


    def log_results(self):
        """
        logs the results of the optimization to an external csv-file.

        :return: None
        """
        # first check, if csv already exists
        if not op.isdir('../Data/Results'):
            os.mkdir('../Data/Results')

        results = {f'current at node {bev}': [self.optimization_model.I[self.current_timestep+1, bev].value]
                   for bev in self.bevs}
        results['timestep'] = [self.current_timestep]

        results = pd.DataFrame(results)

        if not op.isfile('../Data/Results/results.csv'):
            results.to_csv('../Data/Results/results.csv', index=False)

        else:
            # if the file already exist, just append the next incoming results
            results.to_csv('../Data/Results/results.csv', index=False, header=False, mode='a')


    def _store_results(self):
        """
        for usage with rolling horizon: after each optimization, querry
        the results of the first timestep in the current horizon and
        append them to results_SOC/results_I

        :return: None
        """
        for bus in self.bevs:
            # Werte aus model abfragen
            SOC = self.optimization_model.SOC[self.current_timestep, bus].value
            I = self.optimization_model.I[self.current_timestep, bus].value
            # und in Ergebnisliste eintragen
            self.results_SOC[bus].append(SOC)
            self.results_I[bus].append(I)


    def run_optimization_fixed_horizon(self, **kwargs):
        """
        solve the optimization model for the complete considered horizon
        in self.optimzation_model

        :param kwargs: get directly passed to solver
        :return: None
        """
        self.solver_factory.solve(self.optimization_model, tee=kwargs['tee'])
        if type(self)._OPTIONS['log results']:
            self.log_results()


    def run_optimization_rolling_horizon(self, complete_horizon, **kwargs):
        """
        solve multiple optimization models for each mini-horizon in the complete
        considered horizon (and prepare the optimization model for the next
        horizon accordingly)
        1. run the optimization for the current horizon
        2. grab the results from the optimization model
        3. prepare everything for the next timestep

        :param complete_horizon: complete considered horizon [h]
        :param kwargs: get directly passed to solver
        :return: None
        """
        steps = int(complete_horizon * 60/self.resolution)
        for i in range(steps):
            print(i)
            self.run_optimization_fixed_horizon(tee=kwargs['tee'])
            self._store_results()
            self._prepare_next_timestep()
            self._setup_model()


    def plot_grid(self):
        """
        gimmick (and needs networkx): plot the grid

        :return: None
        """
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
        """
        get SOC fullfillment (how much the BEVs have been charged in relation to how
        much they wished to be charged).

        :param optimized: wheather or not a optimization was used
        :return: dict of SOC fullfillments
        """
        final_timestep = self.horizon_width * int(60 / self.resolution) - 1
        if optimized:
            # check if rolling horizon was used
            if not self.rolling:
                # if not, calculate from optimization model Variable instances
                final_socs = {bev: [(self.optimization_model.SOC[final_timestep, bev].value - self.bevs[bev].soc_start) / (self.bevs[bev].soc_target - self.bevs[bev].soc_start) * 100] for bev in self.bevs}
                return final_socs

            else:
                # if yes, calculate from BEVs currents_soc
                final_socs = {bev.home_bus: [(bev.current_soc - bev.soc_start) / (bev.soc_target - bev.soc_start) * 100] for bev in self.bevs.values()}
                return final_socs

        else:
            # if no optimization was used, calculate from BEV current_soc
            final_socs = {bev.home_bus: [(bev.current_soc - bev.soc_start) / (bev.soc_target - bev.soc_start) * 100] for bev in self.bevs.values()}
            return final_socs


    def export_socs_fullfillment(self, optimized=True):
        """
        Export SOC fullfillments as .dat file.

        :param optimized: gets directly passed to _get_soc_fullfillments
        :return: None
        """
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
        """
        for usage in parallel operation with EMO: export results of optimized
        BEV charging currents of the first n timesteps.

        :param width: number of timesteps to consider for export
        :return: dict of optimized BEV charging currents
        """
        ts = self.current_timestep
        return {node: [self.optimization_model.I[time, node].value for time in range(ts, ts+width)] for node in self.bevs}


    @classmethod
    def set_options(cls, key, value):
        """
        override class options to determine how the optimization model
        is going to look (activate aditional constraints or logging of
        data).

        :param key: option name
        :param value: option value
        :return: None
        """
        cls._OPTIONS[key] = value


    def get_results(self, i):
        # function getting called inside the animation. For each call,
        # build and solve the optimization model for the current horizon
        # and return results for first timestep in the horizon
        self.run_optimization_fixed_horizon(tee=False)
        pred = {
            bev: [self.optimization_model.SOC[t, bev].value
                  for t in self.times]
            for bev in self.bevs
        }
        self._store_results()
        self._prepare_next_timestep()
        self._setup_model()
        return {bev: self.results_SOC[bev][i] for bev in self.bevs}, pred


    def animate(self, i, ax, x, ys):
        # first get results (values for first timestep of the current
        # horizon
        res, pred = self.get_results(i)
        x.append(i)
        # y prediction
        #y_pred = self.get_SOC_results()#[self.export_I_results()[bus] for bus in self.bevs]
        print(pred[0])
        x_pred = [i for i in range(self.current_timestep, int(self.current_timestep+self.horizon_width*60/self.resolution))]
        # after the first run (i>0) remove the predictions of the last run
        # (to not get a completely filled plot (in case predictions might
        # change)).
        if i > 0:
            for _ in range(len(self.bevs)):
                ax.lines.pop()

        for num, bev in enumerate(self.bevs):
            ys[num].append(res[bev])
            ax.plot(x, ys[num])

        # add the lines for the prediction of remaining horizon
        for bev in self.bevs:
            ax.plot(x_pred, pred[bev][:], color='gray')

        return ax.lines


    def dyn_gen(self):
        steps = self.horizon_width * 60 / self.resolution
        while self.current_timestep < steps:
            yield self.current_timestep


    def plot_live(self, sudden_load=False):
        # in this function solve subsequentially each horizon and after
        # solution of each horizon plot the data for the fist timestep in the
        # horizon
        fig, ax = plt.subplots(figsize=(15,6))
        ax.set_ylim([0, 110])
        ax.set_xlim([0, (60/self.resolution * self.horizon_width)*2])
        ax.set_xlabel('Zeitschritt')
        ax.set_ylabel('SOC [%]')
        ax.grid()
        ax.plot([96,96], [0, 110])
        x = []
        ys = [[] for _ in self.bevs]
        # a line for course of soc over time for each charger
        #lines = {bev: None for bev in self.bevs}
        #for bev in self.bevs:
            #line, = ax.plot(x, ys)
            #lines[bev.home_bus] = line

        anim = matplotlib.animation.FuncAnimation(
            fig, self.animate, interval=0, fargs=(ax, x, ys), blit=True,
            frames=self.dyn_gen, repeat=False
        )

        plt.show()


    def add_sudden_load(self, start, end, loads_at_buses=None):
        self.sudden_load = SuddenLoad(self, start, end, loads_at_buses)



class SuddenLoad:
    """
    to add sudden loads while running the optimization in rolling horizon.
    Choose in which horizon and at which timesteps there will be which
    additinal load at which bus
    """
    def __init__(
            self,
            glo_object,
            first_horizon,
            last_horizon,
            loads_at_buses=None
    ):
        self.glo_object = glo_object
        self.first_horizon = first_horizon
        self.last_horizon = last_horizon
        self.loads_at_buses = loads_at_buses
        self.prepare_loads_at_buses()
        self.check_sanity()


    def __print__(self):
        return(
            f'Sudden Load effective from timestep {self.first_horizon} '
            f'till timestep {self.last_horizon}.'
        )


    def check_sanity(self):
        glo_buses = self.glo_object.buses
        specified_buses = [bus for bus in self.loads_at_buses]
        if not set(specified_buses).issubset(glo_buses):
            raise ValueError(
                'Specified buses not available in '
                'specified GridLineOptimizer object.'
            )


    def prepare_loads_at_buses(self):
        if self.loads_at_buses != None:
            pass

        else:
            self.loads_at_buses = {
                bus: 0 for bus in self.glo_object.buses
            }


    def set_load_at_all_buses(self, load):
        """
        sets the same load for all present buses

        :param load: load (kW) to be set
        :return: None
        """
        self.loads_at_buses.update(
            {bus: load for bus in self.loads_at_buses}
        )






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

