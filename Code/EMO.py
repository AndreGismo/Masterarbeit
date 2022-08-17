"""
TH Koeln, PROGRESSUS

Module Electrics Management and Optimization (EMO)

Created August 2020

@author: Hotz
00==============================================
Erweitert von André Ulrich für die Masterarbeit
Zusammenspiel von GLO und EMO
What has to be done:
- create grid according to GLO
- pass the household load profiles to GLO
  (then GLO calculates optmized wallbox currents)
- receive the optimized wallbox currents from GLO
"""

#plt.rcParams['figure.dpi'] = 300
import pandas as pd

if True:  # Import
    import pandapower as pp          # grid calculation tool
    import pandapower.plotting       # plotting
    from   scipy.io import loadmat   # load matlab database (household loads)
    import random                    # random numbers
    import numpy                     # vectors, matrices, arrays, gaussians, corrcoeff
    import math                      # various maths stuff
    import cmath                     # complex maths
    import pickle                    # save and load python objects
    import pandas                    # algebra stuff
    import copy
    import sys
    #sys.path.append('C:/Users/Hotz/Dropbox/Weitere/Programmierung')
    #import PSL
    from matplotlib import pyplot as plt
    import matplotlib.dates as mdates
    
    
class Household(): 
    # class household
    def __init__(self,system,position_load,position_bus,position_households,position_main_branch):
        self.position_load = position_load  # position in list "loads" of panda_grid
        self.position_bus = position_bus  # position in list "bus" of panda_grid
        self.position_households = position_households
        self.position_main_branch = position_main_branch
        self.timeline = []
        self.wallbox = None

class Household_Timeline_Handler(): 
    # class for household load data handling, data from IÖW-Verbrauchsdaten
    def __init__(self, system, path, month=1):
        self.system = system
        self.load_data(path)
        self.month = month

    def load_data(self, path): 
        # load and parse household_load_data from MAT-file
        # Pflugradt should be used. 
        tmp = loadmat(path)
        tmp = tmp.get("Verbrauchsdaten").tolist()
        tmp = tmp[0];tmp = tmp[0];tmp = tmp[0]
        tmp = tmp.tolist()  # foo
        tmp = tmp[0];tmp = tmp[0];tmp = tmp[0]
        tmp = tmp.tolist()  # bar
        tmp = tmp[0];tmp = tmp[0];tmp = tmp[0]
        tmp = tmp.tolist()  # tf mate
        tmp = tmp[0]
        self.data = tmp

    def get_random_start_time(self): 
        # gets random day out of specified month to make sure household timelines are different
        day = round(random.random() * 30)
        offset = round(random.random() * 75) # few minutes offset in the time of day for the same reason
        start_frame = ((self.month - 1) * 30 + day) * 1440 + offset
        end_frame = ((self.month - 1) * 30 + day + 1) * 1440 + offset
        return [start_frame, end_frame]

    def make_village_household_loads(self, start_minute, end_minute):  
        # creates timelines and writes them into household instances
        for household_item in self.system.households:
            for index in range(0, int(end_minute / 1430) + 1):  
                # make as many day-timelines as needed plus a litte slack
                data_range = self.get_random_start_time()
                household_item.timeline += self.data[data_range[0] : data_range[1]]
            household_item.timeline = household_item.timeline[start_minute:end_minute]  # crop to required size
            
class Wallbox():  
    # e-vehicle-charger. handles load timelines, SOC, charging. later also controls, communication
    
    WB_index = 0
    
    def __init__(self,system,household,position_load,position_wallboxes,position_bus,
                 rated_power_kW=-1,battery_capacity_kWh=-1,battery_SOC=-1):
        self.index=Wallbox.WB_index
        Wallbox.WB_index+=1
        
        self.rated_power_kW = rated_power_kW
        if (rated_power_kW == -1):  # rated power of charger, if not specified: random value
            self.rated_power_kW = PSL.gaussian_limited(20, 20, 10, 35)
            
        self.battery_capacity_kWh = battery_capacity_kWh  # capacity of battery
        if battery_capacity_kWh == -1:  # if not specified: random value
            self.battery_capacity_kWh = PSL.gaussian_limited(50, 50, 15, 90)  # Tesla Modell S: 100 kWh
        
        self.battery_SOC = battery_SOC
        if battery_SOC == -1:  # initial state of charge; if not specified: random SOC
            self.set_random_SOC()    
        
        self.system             = system  # reference to the system instance it is part of
        self.position_wallboxes = (position_wallboxes)  # position in list "wallboxes" of low_voltage_grid
        self.household          = household  # instance of household it is connected to
        self.household.wallbox  = self
        self.position_load      = position_load  # position in list "loads" of panda_grid
        self.position_bus       = position_bus  # position in list "bus" of panda_grid
        self.ctrl_power_cap     = self.system.max_bus_power  # power cap for voltage controls


    def adapt_ctrl_power_cap(self): 
        # voltage-power-controller. each wallbox has own cap, they're all broadcastet and the global min used
        if min([x for x in self.system.panda_grid.res_bus.vm_pu]) < 0.98:
            self.ctrl_power_cap -= 2
        
        if min([x for x in self.system.panda_grid.res_bus.vm_pu]) > 0.985:
            self.ctrl_power_cap += 0.2
        
        if self.ctrl_power_cap > self.system.max_bus_power:
            self.ctrl_power_cap = self.system.max_bus_power
        
        if self.ctrl_power_cap < 0:
            self.ctrl_power_cap = 0

    def set_random_SOC(self):
        self.battery_SOC = PSL.gaussian_limited(0.3, 0.5, 0, 1)

    def make_random_connection_timeline(self, start_minute, end_minute):  
        # create timeline of when car is connected to charger . to be replaced with work of Sir Marian Sprünken
        
        tmp = []
        for index in range(0, round(end_minute / 1438) + 2):  # number of days the sim runs plus some slack
            # random disconnection time: agent leaves house in the morning    
            disconnect_frame    = int(PSL.gaussian_limited(60 * 7.5, 150, 60 * 4, 60 * 11))  
            # random connection time: agent returns home in the evening
            connect_frame       = int(PSL.gaussian_limited(60 * 17.5, 60 * 4, 60 * 13.5, 60 * 23))  
            tmp += (numpy.ones(disconnect_frame).tolist() 
                    + numpy.zeros(connect_frame - disconnect_frame).tolist() 
                    + numpy.ones(1441 - connect_frame).tolist())
            
        self.connection_timeline = tmp[start_minute:end_minute]  # cut to required length
        
        
    def make_random_connection_timeline_rapid(self, start_minute, end_minute):  
        # create timeline of when car is connected to charger . to be replaced with work of Sir Marian Sprünken
                
        no_boxes=self.system.no_wallboxes
        index=self.index
        
        tmp = []
        tmp += numpy.zeros(index*2+1).tolist()+[1]+numpy.zeros((no_boxes-index-1)*2).tolist()
        loops=math.ceil((end_minute-start_minute)/(len(tmp)))+1
        tmp=tmp*loops
        self.connection_timeline = tmp[0:(end_minute-start_minute+1)]  # cut to required length
        
    def get_power(self, frame, voltage_control):  
        # get power output for frame based on rated power and control algorithm
        
        if (self.connection_timeline[frame] == 1 and self.connection_timeline[frame - 1] == 0): # car arrives at home
                self.set_random_SOC()
        
        if (self.connection_timeline[frame] == 1 and self.battery_SOC < 1):  # car connected and not fully charged
            global_power_cap = min([x.ctrl_power_cap for x in self.system.wallboxes])  # global power cap 
            if voltage_control:
                power = min(self.rated_power_kW, global_power_cap)
            else:
                power = self.rated_power_kW
            self.battery_SOC = (self.battery_SOC + power / 60 / self.battery_capacity_kWh)  # update battery SOC
            return power
        else:
            return 0

class Low_Voltage_System():  
    # class to wrap them all: panda-grid, wallboxes, households
    def __init__(self,line_type,
                 transformer_type):
        self.line_type=line_type
        self.transformer_type=transformer_type
        self.max_bus_power=35
        self.grid = None


    def make_system_from_excel_file(self, file):
        # read stuff from Excel file and format
        lines=pandas.read_excel(file, sheet_name="Lines")
        busses=pandas.read_excel(file, sheet_name="Busses")
        lines=pandas.DataFrame(lines, columns= ['From Bus','To Bus','Length'])
        busses=pandas.DataFrame(busses, columns= ['X','Y','Household','Wallbox'])
        
        household_list     = [pos for (pos,x) in enumerate(busses['Household'],start=0) if x=='Yes']
        wallbox_list       = [pos for (pos,x) in enumerate(busses['Wallbox'],start=0) if x=='Yes']
        
        #wallbox_list= [x for x in sim_handler_1.system.wallboxes if random.random()>0.8]
        
        self.no_households = len(household_list)
        self.no_wallboxes  = len(wallbox_list)
        
        bus_link_list = [list(x) for x in zip(list(lines['From Bus']),list(lines['To Bus']),list(lines['Length']))]
        coordinates   = [list(x) for x in zip(list(busses['X']),list(busses['Y']))]
        
        grid = pp.create_empty_network()
        pp.create_bus(grid, name="ONS_MS", vn_kv=10, type="b")  # bus0 for infinite grid, Medium Voltage
        pp.create_ext_grid(grid, 0, name="starres_Netz", vm_pu=1, va_degree=0)  # define bus as infinite grid
        pp.create_bus(grid, name="ONS_NS", vn_kv=0.4, type="b")  # ONS bus1

        for index in range(2, max([x[1] for x in bus_link_list]) + 1):  # buses on main axes
            pp.create_bus(grid, vn_kv=0.4, name="Bus No. " + str(index))

        for item in bus_link_list: 
            # main grid lines
            pp.create_line(grid,from_bus=item[0],to_bus=item[1],length_km=item[2] * 0.001,
                           std_type=self.line_type,name="Main branch " + str(item))

        self.households = []
        self.wallboxes  = []

        for item in household_list: 
            #create households, wallboxes, connectors etc.
            index_bus_household = pp.create_bus(grid, vn_kv=0.4, name="Bus Household " + str(item))  # bus for house
            pp.create_line(grid,from_bus=item,to_bus=index_bus_household, # line from main axis to house
                           length_km = 10 * 0.001,std_type=self.line_type,
                           name="Main branch " + str(item) + " to household "+ str(index_bus_household))  
            index_load_household = pp.create_load(grid,name="Load household " + str(item),bus=index_bus_household,
                                                  p_mw=0,q_mvar=0)  # create load
            self.households.append(Household(system=self,position_load=index_load_household,
                                             position_bus=index_bus_household,position_households=len(self.households),
                                             position_main_branch=item))

            if item in wallbox_list:
                index_bus_wallbox = pp.create_bus(grid, vn_kv=0.4, name="Bus Wallbox " + str(item))  # bus for WB
                pp.create_line(grid,from_bus=index_bus_household,to_bus=index_bus_wallbox,length_km=1 * 0.001,
                    std_type=self.line_type,name="Household " + str(index_bus_household)+ " to Wallbox "
                    + str(index_bus_wallbox))  # createline from main axis to house
                wallbox_load = pp.create_load(grid,name="Load Wallbox " + str(item),bus=index_bus_wallbox,
                    p_mw=0,q_mvar=0)
                self.wallboxes.append(Wallbox(system=self,household=self.households[-1],position_load=wallbox_load,
                        position_wallboxes=len(self.wallboxes),position_bus=index_bus_wallbox))
                self.households[-1].wallbox = self.wallboxes[-1]
        
        #pp.plotting.create_generic_coordinates(grid)
        pp.runpp(grid)
        
        # for (index,item) in enumerate(coordinates,start=0):
        #     grid.bus_geodata.x[index]=item[1]
        #     grid.bus_geodata.y[index]=item[0]
        #
        # for item in self.households:
        #     grid.bus_geodata.x[item.position_bus]=grid.bus_geodata.x[item.position_main_branch]+0.2
        #     grid.bus_geodata.y[item.position_bus]=grid.bus_geodata.y[item.position_main_branch]
        #
        # for item in self.wallboxes:
        #     grid.bus_geodata.x[item.position_bus]=grid.bus_geodata.x[item.household.position_main_branch]+0.2
        #     grid.bus_geodata.y[item.position_bus]=grid.bus_geodata.y[item.household.position_main_branch]+0.2
        
        pp.create_transformer(grid, 0, 1, name="Transformator", std_type=self.transformer_type)
        return grid


    def grid_from_GLO(self, GLO_grid_file, GLO_grid_params, ideal=True):
        """ By André
        reads in the GLO_grid_file (which contains the information about the
        grid the GLO is optimizing for) and turns it into a pandapower grid
        set as attribute of class
        :param GLO_grid_file: excel-file
        :return: None
        """
        grid_data = pd.read_excel(GLO_grid_file, sheet_name=['Lines', 'Busses'])
        #buses = pd.read_excel(GLO_grid_file, sheet_name='Busses')
        v_mv = 20
        v_lv = 0.4
        s_trafo = GLO_grid_params['S transformer']
        line_specific_impedances = GLO_grid_params['line specific impedances']
        line_lenghts = GLO_grid_params['line lenghts']
        line_capacities = GLO_grid_params['line capacities']
        #i_max_line = GLO_grid_params['i line max']
        # assume all lines to be 15m for calculating R/l
        #r_spec = line_impedance*1000/15 # /1000 since pandas expects length in km

        if ideal:
            vkr = 0
            pfe = 0
            i0 = 0
            vk = 1

        else:
            vkr = 1.5
            pfe = 0.4
            i0 = 0.4
            vk = 6

        # empty grid
        grid = pp.create_empty_network()
        # buses for attaching the transformer
        pp.create_bus(grid, name='transformer mv', vn_kv=v_mv)
        pp.create_bus(grid, name='transformer lv', vn_kv=v_lv)
        # create transformer
        pp.create_transformer_from_parameters(grid, hv_bus=0, lv_bus=1, sn_mva=s_trafo/1000,
                                              vn_hv_kv=v_mv, vn_lv_kv=v_lv, vkr_percent=vkr,
                                              pfe_kw=pfe, i0_percent=i0, vk_percent=vk)
        # slack
        pp.create_ext_grid(grid, bus=0)
        # create all the other busses
        for i in range(2, len(grid_data['Lines'])+2): #+2
            pp.create_bus(grid, name='bus'+str(i), vn_kv=0.4)
            if grid_data['Busses'].loc[i, 'Household'] == 'Yes':
                pp.create_load(grid, bus=i, p_mw=0, name='Load at bus'+str(i))

        # create all the lines
        for i in range(2, len(grid_data['Lines'])+2):
            pp.create_line_from_parameters(grid, from_bus=i-1, to_bus=i, length_km=line_lenghts[i-2]/1000,
                                           r_ohm_per_km=line_specific_impedances[i-2]*1000,
                                           name='line '+str(i-1)+'-'+str(i),
                                           x_ohm_per_km=0, c_nf_per_km=0, max_i_ka=line_capacities[i-2]/1000)

        self.grid = grid



#
class Simulation_Handler():  
    # executes timeseries sim step by step, logs results, handles global voltage-power-controller
    def __init__(self, system, start_minute=0, end_minute=0, rapid=False):
        self.system = system
        self.results = {}    # logged as a dict so it can easily be pickled, to be used e.g. in Topology Estimator
        self.start_minute = start_minute
        self.end_minute = end_minute
        self.rapid=rapid
        # uncomment for data loading from GLO
        #self.create_timelines(self.rapid)  # timelines for households from .mat-file
        # results of the simulation run with GLO (run_GLO_sim) => predicted values
        # starting from the current timestep
        self.res_GLO_sim = {'buses': {i: [] for i in self.system.grid.bus.index if i > 1},
                            'lines': {i: [] for i in self.system.grid.line.index},
                            'trafo': []}

        # store the first values of the results of each timestep here
        # => for plotting afterwards
        self.res_GLO_sim_I = {bus: [] for bus in self.res_GLO_sim['buses']}
        self.res_GLO_sim_U = {line: [] for line in self.res_GLO_sim['lines']}
        self.res_GLO_sim_trafo = []

        
    def create_timelines(self,rapid):
        if rapid:
            self.flat_timelines()
        else:
            self.noisy_timelines()

    def noisy_timelines(self): 
    # timelines with household loads as noise and realistic charging timelines
        for wallbox_item in self.system.wallboxes:
            wallbox_item.make_random_connection_timeline(self.start_minute, self.end_minute)

        self.household_timeline_handler = Household_Timeline_Handler(self.system, household_load_file)
        self.household_timeline_handler.make_village_household_loads(start_minute=self.start_minute, 
                                                                     end_minute=self.end_minute)
    def flat_timelines(self): 
        # noise free timelines for quick and clean CL-matrix generation
        no_of_frames = len(self.system.wallboxes)*2+3 # number of days to be simulated
        for wallbox_item in self.system.wallboxes:
            wallbox_item.make_random_connection_timeline_rapid(self.start_minute, self.end_minute)
    
        self.start_minute = 0
        self.end_minute   = no_of_frames
               
        for item in self.system.households:
            item.timeline = numpy.zeros(no_of_frames ).tolist()
    
    def set_loads(self, voltage_control):  # gets loads from households and wallboxes, writes to PandaGrid
        for household_item in self.system.households:
            self.system.panda_grid.load.p_mw[household_item.position_load] = (
                household_item.timeline[self.frame] * 1e-6)
            self.system.panda_grid.load.q_mvar[household_item.position_load] = (
                household_item.timeline[self.frame] * 0.3 * 1e-6)

        for wallbox_item in self.system.wallboxes:
            self.system.panda_grid.load.p_mw[wallbox_item.position_load] = (
                wallbox_item.get_power(self.frame, voltage_control=voltage_control)* 0.001)
            wallbox_item.adapt_ctrl_power_cap()

    def log_data(self, data, name):  # logs sim results into self.results
        
        if name not in self.results.keys():  # if no item of that name exists yet (first iteration)
            self.results[name] = [] # add new result item to list
        self.results[name].append(data)  # add value to data list of results object, find results item by name

    # hier würden dann die Ergebnisse der Optimierung mit reinkommen => die mal überschreiben => run_sim_stepwise oder sowas
    def run_sim(self, voltage_control=False):  # main loop for timeseries sim
        
        for self.frame in range(1, self.end_minute - self.start_minute - 1):
            #if self.frame / 100 == round(self.frame / 100):
            if self.frame % 100 == 0:
                print("Iteration no. "+ str(self.frame)+ " of "+ str(self.end_minute - self.start_minute))
            
            self.set_loads(voltage_control)
            pp.runpp(self.system.panda_grid)  # LFNR from PandaPower
            
            # gather data that is recorded as simulation results, e.g. bus powers and voltages
            caps    = [x.ctrl_power_cap for x in self.system.wallboxes] # Wallbox power caps from P(U)-controls:
            S       = [x+1j*y for x,y in zip(self.system.panda_grid.load.p_mw, # Power on all busses
                                             self.system.panda_grid.load.q_mvar)] 
            U       = [cmath.rect(self.system.panda_grid.res_bus.vm_pu[x], # Voltages on all busses
                                  self.system.panda_grid.res_bus.va_degree[x]*math.pi/180) 
                       for x in range(len(self.system.panda_grid.res_bus.vm_pu))] 
            u_min   = min([abs(x) for x in U]) # lowest voltage magnitude in the system (for the P(U)-controls)
            I_lines = [1000*sim_handler_1.system.panda_grid.res_line['i_from_ka'][index] # Line current magnitudes
                       for index in range(len(sim_handler_1.system.panda_grid.res_line))]
            
            
            #self.log_data(data= caps,   name="wallbox_caps")  # cap of P(U)-controls
            self.log_data(data= U,       name="U")             # complex voltages on all busses
            self.log_data(data= u_min,   name="u_min")         # complex powers on all busses
            self.log_data(data= S,       name="S")             # complex powers on all busses
            self.log_data(data= I_lines, name="I_lines")       # complex powers on all busses
            
        self.results['pos_U_WB'] = [x.position_bus      # indexes of Wallboxes Voltages in voltages list
                                    for x in self.system.wallboxes]  
        self.results['pos_U_HH'] = [x.position_bus      # indexes of Households Voltages in voltages list 
                                    for x in self.system.households] 
        self.results['pos_S_WB'] = [x.position_load     # indexes of Wallboxes Powers in voltages list 
                                    for x in self.system.wallboxes]  
        self.results['pos_S_HH'] = [x.position_load     # indexes of Households Powers in voltages list 
                                    for x in self.system.households] 
        self.results['Y']        = PSL.get_Y(self.system.panda_grid)
        
        print("Simulation finished (voltage controls="+str(voltage_control)+', rapid='+str(self.rapid)+')')


    def get_results(sim_func):
        def inner(*args, **kwargs):
            sim_func(*args, **kwargs)
            if kwargs['parallel']:
                # args[0] to access self
                args[0].get_first_results()


        return inner


    def get_first_results(self):
        #self.currents = {num: i for num, i in enumerate(self.res_GLO_sim['lines'].values())}
        #self.voltages = {num: u for num, u in enumerate(self.res_GLO_sim['buses'].values())}
        # grab the value of the first timestep of each bus/node and append it
        # to results list (that can late be plotted)
        for line in self.res_GLO_sim['lines']:
            #print(self.res_GLO_sim['lines'])
            #print(line)
            self.res_GLO_sim_I[line+2].append(self.res_GLO_sim['lines'][line][0])

        for bus in self.res_GLO_sim['buses']:
            #print(self.res_GLO_sim['buses'])
            self.res_GLO_sim_U[bus-2].append(self.res_GLO_sim['buses'][bus][0])

        self.res_GLO_sim_trafo.append(self.res_GLO_sim['trafo'][0])



    #@get_results
    def run_GLO_sim(self, household_data, wallbox_data, parallel):#, timesteps, parallel):
        """ By André
        runs a simulation according to the data and grid the GLO was optimizing for
        and stores results (line loading, bus voltage, transformer loading) as an
        attribut of the class, ready to be plotted.
        :param household_data:
        :param wallbox_data:
        :param timesteps:
        :return: None
        """
        self.res_GLO_sim = {'buses': {i: [] for i in self.system.grid.bus.index if i > 1},
                            'lines': {i: [] for i in self.system.grid.line.index},
                            'trafo': []}

        for step in range(len(next(iter(wallbox_data.values())))):
            print('step: ', step)
            # set household loads
            for bus in self.system.grid.load.index:
                self.system.grid.load.loc[bus, 'p_mw'] = household_data[bus][step] * 1e-6

            # add wallbox loads (current*voltage(assumed))
            for wallbox_bus in wallbox_data.keys():
                self.system.grid.load.loc[wallbox_bus, 'p_mw'] += wallbox_data[wallbox_bus][step] *400*1e-6

            # all loads set => start simulation
            pp.runpp(self.system.grid, max_iterations=30)

            # store results
            for line_nr in self.system.grid.line.index:
                self.res_GLO_sim['lines'][line_nr].append(self.system.grid.res_line.loc[line_nr, 'loading_percent'])

            for bus_nr in self.system.grid.bus.index:
                if bus_nr > 1: # skip first two buses
                    self.res_GLO_sim['buses'][bus_nr].append(self.system.grid.res_bus.loc[bus_nr, 'vm_pu']*400)

            self.res_GLO_sim['trafo'].append(self.system.grid.res_trafo.loc[0, 'loading_percent'])

            if parallel:
                # get the results from self.res_GLO_sim[<element>] and append them to
                # self.res_GLO_sim_<element>, so they dont get lost between multiple calls
                self.get_first_results()


    def run_unoptimized_sim(self, household_data, bevs, timesteps, control=False):
        """By André
        run simulation according to unoptimized BEV loading timelines. The generation
        of the loading timeline of each BEV is done inside the BatteryElectricVehicle
        class. These results are passed to the simulation here

        :param household_data: uncontrollable loads from households
        :param bevs: BatteryElectricVehicle objects
        :param timestpes:
        :return:
        """
        cap = 1
        for step in range(timesteps):
            #print(step)
            # set the household loads
            for bus in self.system.grid.load.index:
                self.system.grid.load.loc[bus, 'p_mw'] = household_data[bus][step] * 1e-6
            # query the loading power of each bev at current timestep
            for bev in bevs:
                self.system.grid.load.loc[bev.home_bus, 'p_mw'] += bev.get_current_power(step, cap=cap) * 1e-3

            # now run the simulation
            pp.runpp(self.system.grid, max_iterations=30)

            # store the results of this simulation step
            for line_nr in self.system.grid.line.index:
                self.res_GLO_sim['lines'][line_nr].append(self.system.grid.res_line.loc[line_nr, 'loading_percent'])

            for bus_nr in self.system.grid.bus.index:
                if bus_nr > 1: # skip first two buses
                    self.res_GLO_sim['buses'][bus_nr].append(self.system.grid.res_bus.loc[bus_nr, 'vm_pu']*400)

            self.res_GLO_sim['trafo'].append(self.system.grid.res_trafo.loc[0, 'loading_percent'])

            # call the controller P(U)
            if control:
                #cap = self.control_power(thres=0.98, cap=0.5)
                cap = self.alt_control_power(thres_lo=0.94, thres_hi=0.98, last_cap=cap)


    def control_power(self, thres, cap):
        """
        run P(U) controlling. Check the voltages at each node, if one falls below critical value
        than lower each controllable load by a certain amount
        :param thres: critical value for node voltage
        :param cap: amount to cap loads
        :return:
        """
        # check all the node voltages
        failed = False
        for bus in self.system.grid.res_bus.index:
            if self.system.grid.res_bus.loc[bus, 'vm_pu'] < thres:
                #print(f"voltage too low: {self.system.grid.res_bus.loc[bus, 'vm_pu']*400}")
                failed = True
                #break

        if failed:
            return 1-cap
        else:
            return 1


    def alt_control_power(self, thres_lo, thres_hi, last_cap):
        if min([x for x in self.system.grid.res_bus.vm_pu]) < thres_lo:
            new_cap = last_cap * (1 - 0.05)
            return new_cap

        elif min([x for x in self.system.grid.res_bus.vm_pu]) > thres_hi:
            new_cap = last_cap * (1 + 0.02)
            if new_cap >= 1:
                new_cap = 1
            return new_cap

        else:
            return last_cap


    def plot_EMO_sim_results(self, freq, element='buses', legend=True, marker='o', save=False,
                             usetex=False, compact_x=False, name_ext='', compare=False):
        if usetex:
            plt.rcParams['text.usetex'] = True
            plt.rcParams['font.family'] = 'serif'
            plt.rcParams['grid.linewidth'] = 0.4
            plt.rcParams['lines.linewidth'] = 1
            plt.rcParams['legend.fontsize'] = 8
            plt.rcParams['font.size'] = 11

        if compact_x:
            x_fmt = mdates.DateFormatter('%H')

        fig, ax = plt.subplots(1, 1, figsize=(6.5, 1.75))
        length = len(self.res_GLO_sim['trafo'])
        date_range = pd.date_range(start='2022', freq=str(freq)+'min', periods=length)
        if element == 'buses':
            for elm_nr in self.res_GLO_sim[element]:
                ax.plot(date_range, self.res_GLO_sim[element][elm_nr], label='Spannung an Bus {}'.format(elm_nr-1),
                         marker=marker)
            ax.set_ylabel('Knotenspannung [V]')
            if compact_x:
                ax.xaxis.set_major_formatter(x_fmt)
                ax.set_xlabel('Time [hh]')
            else:
                ax.set_xlabel('Zeitpunkt [MM-TT hh]')
            #ax.set_title('Knotenspannungen über der Zeit', fontsize=20)
            ax.grid()
            if legend:
                ax.legend()

            if save:
                fig.savefig(f'res_sim_buses{name_ext}.pdf', bbox_inches='tight')

        elif element == 'lines':
            for elm_nr in self.res_GLO_sim[element]:
                ax.plot(date_range, self.res_GLO_sim[element][elm_nr], label='Auslastung von Leitung {}'.format(elm_nr),
                         marker=marker)
            if usetex:
                ax.set_ylabel(r'Auslastung [\%]')

            else:
                ax.set_ylabel('Auslastung [%]')

            if compact_x:
                ax.xaxis.set_major_formatter(x_fmt)
                ax.set_xlabel('Time [hh]')
            else:
                ax.set_xlabel('Zeitpunkt [MM-TT hh]')
            #ax.set_title('Auslastung der Leitungen über der Zeit', fontsize=20)
            ax.grid()
            if legend:
                ax.legend()

            if save:
                fig.savefig(f'res_sim_lines{name_ext}.pdf', bbox_inches='tight')

        elif element == 'trafo':
            if not compare:
                ax.plot(date_range, self.res_GLO_sim[element], marker=marker)
            else:
                ax.plot(date_range, self.res_GLO_sim[element], marker=marker, label='ohne Optimierung')
                ax.plot(date_range, self.last_results[element], marker=marker, label='mit Optimierung')
                if legend:
                    ax.legend()

            if usetex:
                ax.set_ylabel(r'Transformer load [\%]')

            else:
                ax.set_ylabel('Auslastung [%]')
            if compact_x:
                ax.xaxis.set_major_formatter(x_fmt)
                ax.set_xlabel('Time [hh]')
            else:
                ax.set_xlabel('Zeitpunkt [MM-TT hh]')
            #ax.set_title('Auslastung des Transformators über der Zeit', fontsize=20)
            ax.grid()

            if save:
                fig.savefig(f'res_sim_trafo{name_ext}.pdf', bbox_inches='tight')


    def export_sim_results(self, element, len_hrs=24, res_min=15, sep='\t'):
        """
        export values of simulated elements for nice plotting with pgfplots
        :param element: which elements values to export
        :param sep: column separator
        :return: none
        """
        data_df = pd.DataFrame(self.res_GLO_sim[element])
        data_df.index = numpy.linspace(0, len_hrs, len_hrs * int(60 / res_min))
        data_df.to_csv(f'export-sim-results-{element}.dat', sep=sep)


    def store_sim_results(self, name_extension=''):
        """
        store the results of the previous simulation as csv.
        :return: None
        """
        buses_df = pd.DataFrame(self.res_GLO_sim['buses'])
        lines_df = pd.DataFrame(self.res_GLO_sim['lines'])
        trafo_df = pd.DataFrame(self.res_GLO_sim['trafo'])

        # save the to csv
        buses_df.to_csv(f'results/res_buses{name_extension}.csv', index=False)
        lines_df.to_csv(f'results/res_lines{name_extension}.csv', index=False)
        trafo_df.to_csv(f'results/res_trafo{name_extension}.csv', index=False)


    def load_sim_results(self):
        buses_df = pd.read_csv('res_buses.csv')
        lines_df = pd.read_csv('res_lines.csv')
        trafo_df = pd.read_csv('res_trafo.csv')

        # store as dict, similar to self.res_GLO_sim
        self.last_results = {'buses': {int(i): [buses_df.loc[j, i] for j in buses_df.index] for i in buses_df.columns},
                             'lines': {int(i): [lines_df.loc[j, i] for j in lines_df.index] for i in lines_df.columns},
                             'trafo': [trafo_df.loc[j, '0'] for j in trafo_df.index]}


    def reset_GLO_sim_results(self):
        # reset the res_GLO_sim dict
        self.res_GLO_sim = {'buses': {i: [] for i in self.system.grid.bus.index if i > 1},
                            'lines': {i: [] for i in self.system.grid.line.index},
                            'trafo': []}



    def plot_powers(self):
        for bus_no in range(len(sim_handler_1.results['S'][0])):
            PSL.plot_list([sim_handler_1.results['S'][frame][bus_no] for frame in range(len(sim_handler_1.results['S']))],
                  title='timeline power bus no. '+str(bus_no))
        
# MAIN-Routine==================================================================================================
if __name__ == '__main__':
    #grid_name='Example_Grid'
    grid_name='Example_Grid_Simple'
    #grid_name='selfmade_grid'

    #grid_name='Example_Grid_Current_Estimation'

    excel_file='grids/'+grid_name+'.xlsx'
    save_file='sav/'+grid_name+'.pic'
    results_dict_file='sav/'+grid_name+'_results.pic'
    household_load_file="zeitreihen/Verbrauchsdaten_IOW"
    handler_file='sav/'+grid_name+'_handler.pic'
    panda_file='sav/'+grid_name+'_panda.pic'
    rapid=False

    if rapid==True:
        results_dict_file='sav/'+grid_name+'_results_rapid.pic'

    if True:  # Simulation rather than loading from file
        system_1 = Low_Voltage_System(line_type='NAYY 4x120 SE',transformer_type="0.25 MVA 10/0.4 kV")
        system_1.panda_grid = system_1.make_system_from_excel_file(file=excel_file)
        system_1.grid_from_GLO('grids/selfmade_grid.xlsx', {'S transformer': 240,
                                                            'line impedance': 0.004,
                                                            'i line max': 140
                                                            })
        # simulate uncontrolled system:
        sim_handler_1 = Simulation_Handler(system_1,
                                           start_minute=60 * 12,
                                           end_minute=60 * 12 + 24 * 60,
                                           rapid=rapid)

        if False: # Simulation w/ voltage controls
            sim_handler_2 = copy.deepcopy(sim_handler_1) # run voltage controlled system simulation
            sim_handler_2.run_sim(voltage_control=True)  # run voltage controlled system simulation
        sim_handler_1.run_sim(voltage_control=False)

        if True: # Save simulation results
            with open(results_dict_file, "wb") as f:  # save to file
                pickle.dump(sim_handler_1.results, f)     # write to file



            with open(panda_file, "wb") as f:  # save to file
                pickle.dump(system_1.panda_grid, f)     # write to file
            if True:
                with open(handler_file, "wb") as f:  # save to file
                    #pickle.dump(sim_handler_1, f)     # write to file
                    pass # uncomment above to avoid error

    else:  # Load simulation results from file rather than run simulation
        with open(results_dict_file,    "rb")   as f:
            res = pickle.load(f)    # read from file
        with open(panda_file,           "rb")   as f:  # save to file
            pg_orig = pickle.load(f)     # write to file
        if True:
            with open(handler_file,         "rb")   as f:
                sim_handler_1 = pickle.load(f)    # read from file

    sim_handler_1.plot_powers()
    # x=[sim_handler_1.results['u_min'], sim_handler_2.results['u_min']]
    # import numpy as np
    # x=np.array(x).T.tolist()
    # PSL.plot_list(x)

    #[1000*sim_handler_1.system.panda_grid.res_line['i_from_ka'][index] for index in range(len(sim_handler_1.system.panda_grid.res_line))]


    # S_WB=sim_handler_1.results['pos_S_WB']
    # S_HH=sim_handler_1.results['pos_S_HH']

    # PSL.plot_list([sum([frame[WB]*1000 for WB in sim_handler_1.results['pos_S_WB']]) for frame in sim_handler_1.results['S']],title='E-Mobilitäts-Lastprofil, kumuliert, 30 Fahrzeuge')
    # PSL.plot_list([sum([frame[WB]*1000 for WB in sim_handler_2.results['pos_S_WB']]) for frame in sim_handler_2.results['S']],title='E-Mobilitäts-Lastprofil, kumuliert, 30 Fahrzeuge, geregelt')
    # PSL.plot_list([[frame[HH] for HH in sim_handler_1.results['pos_S_HH']] for frame in sim_handler_1.results['S']])
    # PSL.plot_list([[frame[HH]*1000 for HH in [6,8]] for frame in sim_handler_1.results['S']],title="Haushaltslast")















