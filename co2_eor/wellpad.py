#wellpad model V2

import pyomo.environ as pyo
import pandas as pd
from pyomo.common.config import ConfigBlock, ConfigValue, In
from idaes.core import (
    ControlVolume0DBlock,
    declare_process_block_class,
    UnitModelBlockData,
    useDefault,
)
from idaes.core.util.config import is_physical_parameter_block
from idaes.core.util.scaling import set_scaling_factor
from idaes.core.scaling.autoscaling import AutoScaler

#specify configuration options
def make_wellpad_config_block(config):
    config.declare(
            "property_package",
            ConfigValue(default=useDefault, domain=is_physical_parameter_block),
            )
    config.declare(
            "has_phase_equilibrium",
            ConfigValue(default=False, domain=In([False]))
            )
    config.declare(
            "temperature",
            ConfigValue(default=300,domain=float) #reservoir temperature
            )
    config.declare(
            "pressure",
            ConfigValue(default=120*100000,domain=float) #reservoir pressure
            )
    config.declare(
            "Boi",
            ConfigValue(default=1.3,domain=float)#RBoil/STBoil
            )
    config.declare(
            "GOR",
            ConfigValue(default=1.5,domain=float)#RBgas/STBoil
            )
    config.declare(
            "PC_A",
            ConfigValue(default=0.4,domain=float)#production curve fit parameter A
            )
    config.declare(
            "PC_B",
            ConfigValue(default=0.7,domain=float)#production curve fit parameter B
            )
    config.declare(
            "GB_A",
            ConfigValue(default=0.6,domain=float)#breakthrough curve fit parameter A
            )
    config.declare(
            "GB_B",
            ConfigValue(default=1.13,domain=float)#breakthrough curve fit parameter B
            )
    config.declare(
            "kovr",
            ConfigValue(default=1.1e-13,domain=float)#overall permeability constant STB
            )
    config.declare(
            "use_correction_factor",
            ConfigValue(default=False,domain=bool)
            )
    config.declare(
            "SC_A",
            ConfigValue(default=0.4,domain=float)#sensitivity curve fit parameter A
            )
    config.declare(
            "SC_B",
            ConfigValue(default=789.2,domain=float)#sensitivity curve fit parameter B
            )
    config.declare(
            "IR_base",
            ConfigValue(default=0.00231,domain=float)#base case injection rate #RB/s
            )
    config.declare(
            "BHP_max",
            ConfigValue(default=1000*100000,domain=float)#maximum injection pressure
            )
    config.declare(
            "multiplier",
            ConfigValue(default=1,domain=int) #flow multiplier
            )
    config.declare(
            "depth",
            ConfigValue(default=1000,domain=float) #reservoir depth
            )
    config.declare(
            "property_package_args",
            ConfigBlock(implicit=True)
            )

#create control volume and required state blocks
def make_control_volume(unit,name,config):
    control_volume = ControlVolume0DBlock(
            property_package=config.property_package,
            property_package_args=config.property_package_args,
            )

    setattr(unit,name,control_volume)
    control_volume.add_state_blocks(has_phase_equilibrium=config.has_phase_equilibrium)
    #create injection state block - if want to use time indexing this will need more thought
    control_volume.injection_state=config.property_package.build_state_block(defined_state=False)
    control_volume.reservoir_state=config.property_package.build_state_block(defined_state=False)

#adding parameters from config
def add_params(unit,name,config):
    unit.reservoir_temperature = pyo.Param(initialize=config.temperature)
    unit.reservoir_pressure = pyo.Param(initialize=config.pressure)
    unit.Boi = pyo.Param(initialize=config.Boi)
    unit.GOR = pyo.Param(initialize=config.GOR)
    unit.PC_A = pyo.Param(initialize=config.PC_A)
    unit.PC_B = pyo.Param(initialize=config.PC_B)
    unit.GB_A = pyo.Param(initialize=config.GB_A)
    unit.GB_B = pyo.Param(initialize=config.GB_B)
    unit.kovr = pyo.Param(initialize=config.kovr)
    if config.use_correction_factor:
        unit.SC_A = pyo.Param(initialize=config.SC_A)
        unit.SC_B = pyo.Param(initialize=config.SC_B)
        unit.IR_base = pyo.Param(initialize=config.IR_base)
        unit.SC_at_base_IR = pyo.Param(initialize=unit.SC_A*(1-pyo.exp(-unit.SC_B*unit.IR_base)))
    unit.BHP_max = pyo.Param(initialize=config.BHP_max)
    unit.multiplier = pyo.Param(initialize=config.multiplier)

#adding variables and constraints
def add_equations(unit,name,config):
    inlet=unit.control_volume.properties_in[0]
    injection_state=unit.control_volume.injection_state
    reservoir_state=unit.control_volume.reservoir_state
    outlet=unit.control_volume.properties_out[0]

    #fix reservoir temperature and pressure
    unit.reservoir_temperature_constraint = pyo.Constraint(
        expr=reservoir_state.temperature==unit.reservoir_temperature
    )
    unit.reservoir_pressure_constraint = pyo.Constraint(
        expr=reservoir_state.pressure==unit.reservoir_pressure
    )

    #define HCPV
    unit.HCPV = pyo.Var(domain=pyo.NonNegativeReals)

    #create slack variables for initialization
    num_slacks=2
    unit.spos = pyo.Var(range(1,num_slacks+1),domain=pyo.NonNegativeReals,initialize=0)
    unit.sneg = pyo.Var(range(1,num_slacks+1),domain=pyo.NonNegativeReals,initialize=0)
    unit.spos.fix(0)
    unit.sneg.fix(0)
    #create feasibility expression and objective function
    unit.feasibility_expression = pyo.Expression(expr=pyo.quicksum(unit.spos[i]+unit.sneg[i] for i in range(1,num_slacks+1)))
    unit.feasibility_objective = pyo.Objective(expr=unit.feasibility_expression)
    unit.feasibility_objective.deactivate()

    #equality constraints

    #connect inlet state and injection state
    unit.inlet_injection_mass_bal = pyo.Constraint(
        expr=inlet.flow_mass==injection_state.flow_mass*unit.multiplier+unit.spos[1]-unit.sneg[1]
    )
    unit.inlet_injection_temperature_bal = pyo.Constraint(
        expr=inlet.temperature==injection_state.temperature
    )
    unit.inlet_injection_pressure_bal = pyo.Constraint(
        expr=inlet.pressure==injection_state.pressure
    )

    #mass balance between injection state and reservoir state
    unit.injection_reservoir_mass_bal = pyo.Constraint(
            expr=injection_state.flow_mass==reservoir_state.flow_mass
            ) #constraint 1

    """
    #CO2 Injection Rate RB/s
    unit.q_CO2_INJ = pyo.Var(domain=pyo.NonNegativeReals)
    #calculate injection rate
    unit.injection_rate_from_density_constraint = pyo.Constraint(
            expr=unit.q_CO2_INJ==injection_state.flow_mass/injection_state.dens_mass*6.29 +unit.spos[1]-unit.sneg[1]
            ) #constraint 3
    unit.injection_rate_from_darcys_law_constraint = pyo.Constraint(
            expr=unit.q_CO2_INJ==
                unit.kovr*6.29/reservoir_state.visc_d_phase["Liq"]*(injection_state.pressure-reservoir_state.pressure)+
                    unit.spos[2]-unit.sneg[2]
            ) #constraint 4
    """
    #combine above to say that injection state flow_vol (m3/s) == dacrys law eqn (m3/s)
    unit.darcys_law = pyo.Constraint(
        expr=injection_state.flow_vol==
                unit.kovr/reservoir_state.visc_d_phase["Liq"]*(injection_state.pressure-reservoir_state.pressure)+
                        unit.spos[2]-unit.sneg[2]
    )

    #sensitivity curve correction factor
    if config.use_correction_factor:
        unit.correction_factor = pyo.Expression(
                expr=(1-pyo.exp(-unit.SC_B*reservoir_state.flow_vol*6.29))/(1-pyo.exp(-unit.SC_B*unit.IR_base))
                )
    else:
        unit.correction_factor = pyo.Expression(expr=1)
    
    #expressions for production rates
    #derivative of production curve, change in incremental recovery factor over change in HCPV
    unit.dRfdHCPV = pyo.Expression(
        expr=unit.PC_A*unit.PC_B/(unit.HCPV+unit.PC_B)**2
    )

    #Oil production rate STB/s
    unit.q_OIL_PROD = pyo.Expression(
            expr=1/unit.Boi*unit.dRfdHCPV*reservoir_state.flow_vol*6.29*unit.correction_factor*unit.multiplier
            )
    #gas production rate
    unit.q_GAS_PROD = pyo.Expression(
            expr=unit.GOR*unit.q_OIL_PROD*unit.Boi*unit.multiplier
            )
    #derivative of gas breakthrough curve
    unit.dGbdHCPV = pyo.Expression(
        expr=unit.GB_A*unit.GB_B*unit.HCPV**(unit.GB_B-1)
    )
    #gas breakthrough rate
    unit.q_CO2_BRKTH = pyo.Expression(
            expr=unit.dGbdHCPV*reservoir_state.flow_vol*6.29*unit.multiplier
            )
    
    #connect reservoir state to outlet state
    unit.reservoir_outlet_mass_bal = pyo.Constraint(
        expr=reservoir_state.flow_mass*unit.multiplier*unit.dGbdHCPV==outlet.flow_mass
    )
    unit.reservoir_outlet_temperature_bal = pyo.Constraint(
        expr=reservoir_state.temperature==outlet.temperature
    )
    unit.reservoir_outlet_pressure_bal = pyo.Constraint(
        expr=reservoir_state.pressure==outlet.pressure
    )
    
    #inequality constraints
    #maximum injection pressure
    unit.max_BHP_constraint = pyo.Constraint(
            expr=injection_state.pressure<=unit.BHP_max
            )
    #ensure injection pressure is greater than reservoir pressure
    unit.min_BHP_constraint = pyo.Constraint(
            expr=injection_state.pressure>=reservoir_state.pressure
            )

def guess_scales(unit,name,config):
    inlet=unit.control_volume.properties_in[0]
    injection_state=unit.control_volume.injection_state
    reservoir_state=unit.control_volume.reservoir_state
    outlet=unit.control_volume.properties_out[0]

    #parameter scales
    set_scaling_factor(unit.reservoir_temperature,1e-2)
    set_scaling_factor(unit.reservoir_pressure,1e-7)
    set_scaling_factor(unit.PC_A,10)
    set_scaling_factor(unit.PC_B,10)
    set_scaling_factor(unit.GB_A,10)
    set_scaling_factor(unit.kovr,1e13)
    if config.use_correction_factor:
        set_scaling_factor(unit.SC_A,10)
        set_scaling_factor(unit.SC_B,1e-2)
        set_scaling_factor(unit.IR_base,1e3)
    set_scaling_factor(unit.BHP_max,1e-7)
    
    #variable scales
    set_scaling_factor(inlet.pressure,1e-7)
    set_scaling_factor(reservoir_state.pressure,1e-7)
    set_scaling_factor(injection_state.pressure,1e-7)
    set_scaling_factor(outlet.pressure,1e-7)
    set_scaling_factor(inlet.temperature,1e-2)
    set_scaling_factor(reservoir_state.temperature,1e-2)
    set_scaling_factor(injection_state.temperature,1e-2)
    set_scaling_factor(outlet.temperature,1e-2)
    set_scaling_factor(inlet.flow_mass,1e1)
    set_scaling_factor(reservoir_state.flow_mass,1e1)
    set_scaling_factor(injection_state.flow_mass,1e1)
    set_scaling_factor(outlet.flow_mass,1e1)

    #equality constraint scales
    set_scaling_factor(unit.reservoir_temperature_constraint,1e-2)
    set_scaling_factor(unit.reservoir_pressure_constraint,1e-7)
    set_scaling_factor(unit.inlet_injection_mass_bal,1e1)
    set_scaling_factor(unit.inlet_injection_temperature_bal,1e-2)
    set_scaling_factor(unit.inlet_injection_pressure_bal,1e-7)
    set_scaling_factor(unit.injection_reservoir_mass_bal,1e1)
    set_scaling_factor(unit.reservoir_outlet_mass_bal,1e1)
    set_scaling_factor(unit.reservoir_outlet_temperature_bal,1e-2)
    set_scaling_factor(unit.reservoir_outlet_pressure_bal,1e-7)
    #set_scaling_factor(unit.injection_rate_from_density_constraint,1e4)
    #set_scaling_factor(unit.injection_rate_from_darcys_law_constraint,1e4)
    set_scaling_factor(unit.darcys_law,1e4)

    #inequality constraint scales
    set_scaling_factor(unit.min_BHP_constraint,1e-7)
    set_scaling_factor(unit.max_BHP_constraint,1e-7)

#define wellpad class
@declare_process_block_class("wellpad")
class wellpadData(UnitModelBlockData):
    CONFIG = UnitModelBlockData.CONFIG()
    make_wellpad_config_block(CONFIG)

    def build(self):
        super(wellpadData,self).build()
        make_control_volume(self,"control_volume",self.config)
        add_params(self,"params",self.config)
        add_equations(self,"constraints",self.config)
        guess_scales(self,'scales',self.config)
        self.add_inlet_port(block=self.control_volume.properties_in,name="inlet")
        self.add_outlet_port(block=self.control_volume.properties_out,name="outlet")
    
    def activate_slack_variables(self):
        self.spos.unfix()
        self.sneg.unfix()

    def deactivate_slack_variables(self):
        self.spos.fix(0)
        self.sneg.fix(0)

    def activate_feasibility_problem(self):
        self.activate_slack_variables()
        self.feasibility_objective.activate()
    
    def deactivate_feasibility_problem(self):
        self.deactivate_slack_variables()
        self.feasibility_objective.deactivate()

    def initialize(self,solver=None,tee=False,display_after=False):
        print(f'initializing {self.name}')

        #activate feasibility problem
        self.activate_feasibility_problem()
        #scale model
        scaled_self = pyo.TransformationFactory('core.scale_model').create_using(self)
        if solver==None:
            solver = pyo.SolverFactory('ipopt')
            solver.options['linear_solver']='ma97'
        res = solver.solve(scaled_self,tee=tee)
        #undo scaling
        pyo.TransformationFactory('core.scale_model').propagate_solution(scaled_self,self)
        if display_after: 
            self.display()
            self.print_all()
        #deactivate feasibility problem
        self.deactivate_feasibility_problem()
        #create autoscaler
        autoScaler=AutoScaler(overwrite=True)
        autoScaler.scale_variables_by_magnitude(self)
        #autoScaler.scale_constraints_by_jacobian_norm(self)

        print(f'{self.name} initialization complete')

    def print_parameters(self):
        print(f'{self.name} parameters')
        print(f'oil volume factor: {pyo.value(self.Boi)}')
        print(f'gas oil ratio: {pyo.value(self.GOR)}')
        print(f'production curve fit parameter A: {pyo.value(self.PC_A)}')
        print(f'production curve fit parameter B: {pyo.value(self.PC_B)}')
        print(f'gas breakthrough curve fit parameter A: {pyo.value(self.GB_A)}')
        print(f'gas breakthrough curve fit parameter B: {pyo.value(self.GB_B)}')
        print(f'kovr: {pyo.value(self.kovr)}')
        print(f'max BHP: {pyo.value(self.BHP_max)}')
        print(f'multiplier: {pyo.value(self.multiplier)}')
        print()

    def print_variables(self):
        print(f'{self.name} variables')
        print(f'HCPV: {pyo.value(self.HCPV)}')
        print(f'CO2 injection rate (kg/s): {pyo.value(self.inlet.flow_mass[0])}')
        print(f'Inlet pressure: {pyo.value(self.inlet.pressure[0])/100000}')
        print(f'Inlet temperature: {pyo.value(self.control_volume.properties_in[0].temperature)}')
        print(f'Injection pressure: {pyo.value(self.control_volume.injection_state.pressure)/100000}')
        print(f'Injection temperature: {pyo.value(self.control_volume.injection_state.temperature)}')
        print(f'Reservoir pressure: {pyo.value(self.control_volume.reservoir_state.pressure)/100000}')
        print(f'Reservoir temperature: {pyo.value(self.control_volume.reservoir_state.temperature)}')
        print()


    def print_expressions(self):
        print(f'{self.name} expressions')
        print(f'correction factor: {pyo.value(self.correction_factor)}')
        print(f'oil production rate: {pyo.value(self.q_OIL_PROD)}')
        print(f'gas production rate: {pyo.value(self.q_GAS_PROD)}')
        print(f'CO2 breakthrough rate: {pyo.value(self.q_CO2_BRKTH)}')
        print()

    def print_all(self):
        print()
        self.print_parameters()
        self.print_variables()
        self.print_expressions()

    def export_df(self):
        data = {
                "reservoir temperature (K)":pyo.value(self.reservoir_temperature),
                "reservoir pressure (bar)": pyo.value(self.reservoir_pressure)/100000,
                "oil volume factor (STB oil/RB oil)":pyo.value(self.Boi),
                "gas oil ratio (RB gas/RB oil)":pyo.value(self.GOR),
                "production curve fit parameter A":pyo.value(self.PC_A),
                "production curve fit parameter B":pyo.value(self.PC_B),
                "gas breakthrough curve fit parameter A":pyo.value(self.GB_A),
                "gas breakthrough curve fit parameter B":pyo.value(self.GB_B),
                "k overall (m3)":pyo.value(self.kovr),
                "max BHP (bar)":pyo.value(self.BHP_max)/100000,
                "multiplier":pyo.value(self.multiplier),
                "HCPV":pyo.value(self.HCPV),
                "total CO2 injection rate (kg/s)":pyo.value(self.inlet.flow_mass[0]),
                "unit CO2 injection rate (kg/s)":pyo.value(self.control_volume.injection_state.flow_mass),
                "unit CO2 injection rate (Rb/s)":pyo.value(self.control_volume.injection_state.flow_vol*6.29),
                "inlet pressure (bar)":pyo.value(self.inlet.pressure[0])/100000,
                "inlet temperature (K)":pyo.value(self.control_volume.properties_in[0].temperature),
                "injection pressure (bar)":pyo.value(self.control_volume.injection_state.pressure)/100000,
                "injection temperature (K)":pyo.value(self.control_volume.injection_state.temperature),
                "correction factor":pyo.value(self.correction_factor),
                "production curve slope":pyo.value(self.dRfdHCPV),
                "gas breakthrough curve slope":pyo.value(self.dGbdHCPV),
                "oil production rate (STB/s)":pyo.value(self.q_OIL_PROD),
                "gas production rate (RB/s)":pyo.value(self.q_GAS_PROD),
                "CO2 breakthrough rate (RB/s)":pyo.value(self.q_CO2_BRKTH),
                "total CO2 out (kg/s)":pyo.value(self.outlet.flow_mass[0])
                }

        return pd.DataFrame(data,index=[self.name])

