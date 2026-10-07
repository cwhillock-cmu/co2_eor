#well pattern model V3

import pyomo.environ as pyo
import pandas as pd
from pyomo.common.config import ConfigBlock, ConfigValue, In
import idaes
from idaes.core import (
    ControlVolume0DBlock,
    declare_process_block_class,
    UnitModelBlockData,
    useDefault,
)
from idaes.core.util.config import is_physical_parameter_block
from idaes.core.util.scaling import set_scaling_factor
from idaes.core.scaling.autoscaling import AutoScaler
from idaes.core.util.math import smooth_max
from idaes.core.util.initialization import propagate_state
from co2_eor.util_funcs import add_state_material_balances
"""
well pattern injection model
pure co2 injection
output is a gas mixture of co2, and currently only ch4
unforunately need to hardcode in the assumption about what equation of states are being used
inlet port uses the helmholtz equation of state 
outlet port uses the modular properties framework
both eos with amount basis mass
"""

#specify configuration options
def make_wellpattern_config_block(config):
    config.declare(
            "primary_property_package",
            ConfigValue(default=useDefault, domain=is_physical_parameter_block),
            )
    config.declare(
            "secondary_property_package",
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
            ConfigValue(default=1.5,domain=float) #Sm3 gas/STBoil
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
            "injectivity_index",
            ConfigValue(default=1.1e-9,domain=float) #injectivity index - another layer of abstraction
    )
    config.declare(
            "reference_pressure",
            ConfigValue(default=120*100000,domain=float) #reference pressure for injectivity equation - abstracted
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
            "pressure_max",
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
            "outlet_pressure",
            ConfigValue(default=40*100000,domain=float) #gasses outlet pressure
    )
    config.declare(
            "outlet_temperature",
            ConfigValue(default=300,domain=float) #gasses outlet temperature
    )
    config.declare(
            "primary_property_package_args",
            ConfigBlock(implicit=True)
            )
    config.declare(
            "secondary_property_package_args",
            ConfigBlock(implicit=True)
            )

#required state blocks
def make_control_volume(unit,name,config):
    #not using flowsheet time for now
    unit.control_volume = pyo.Block()

    #create inlet state block
    #uses primary property package, transitioning to cubic eos
    unit.control_volume.properties_in = config.primary_property_package.build_state_block(unit.flowsheet().time,defined_state=True,has_phase_equilibrium=config.has_phase_equilibrium)

    #reservoir state
    unit.control_volume.reservoir_state = config.primary_property_package.build_state_block(unit.flowsheet().time,defined_state=True,has_phase_equilibrium=config.has_phase_equilibrium)

    #create outlet state block
    #uses secondary property package, mixture of gasses
    unit.control_volume.properties_out = config.secondary_property_package.build_state_block(unit.flowsheet().time,defined_state=False,has_phase_equilibrium=config.has_phase_equilibrium)

#adding parameters from config
def add_params(unit,name,config):
    unit.reservoir_temperature = pyo.Param(initialize=config.temperature) #K
    unit.reservoir_pressure = pyo.Param(initialize=config.pressure) #Pa
    set_scaling_factor(unit.reservoir_temperature,1e-2)
    set_scaling_factor(unit.reservoir_pressure,1e-7)

    unit.Boi = pyo.Param(initialize=config.Boi) #RB oil / STB oil
    unit.GOR = pyo.Param(initialize=config.GOR) #RB gas / STB oil

    unit.PC_A = pyo.Param(initialize=config.PC_A) 
    unit.PC_B = pyo.Param(initialize=config.PC_B)
    unit.GB_A = pyo.Param(initialize=config.GB_A)
    unit.GB_B = pyo.Param(initialize=config.GB_B)
    set_scaling_factor(unit.PC_A,10)
    set_scaling_factor(unit.PC_B,10)
    set_scaling_factor(unit.GB_A,10)
    
    unit.kovr = pyo.Param(initialize=config.kovr) #m^3
    set_scaling_factor(unit.kovr,1e13)
    unit.injectivity_index = pyo.Param(initialize=config.injectivity_index) #m^3/Pa
    set_scaling_factor(unit.injectivity_index,1e9)
    unit.reference_pressure = pyo.Param(initialize=config.reference_pressure) #Pa
    set_scaling_factor(unit.reference_pressure,1e-7)

    if config.use_correction_factor:
        unit.SC_A = pyo.Param(initialize=config.SC_A)
        unit.SC_B = pyo.Param(initialize=config.SC_B)
        unit.IR_base = pyo.Param(initialize=config.IR_base) #RB/s -- double check
        unit.SC_at_base_IR = pyo.Param(initialize=unit.SC_A*(1-pyo.exp(-unit.SC_B*unit.IR_base)))
        set_scaling_factor(unit.SC_A,10)
        set_scaling_factor(unit.SC_B,1e-2)
        set_scaling_factor(unit.IR_base,1e3)
    
    unit.pressure_max = pyo.Param(initialize=config.pressure_max) #Pa
    set_scaling_factor(unit.pressure_max,1e-7)
    unit.multiplier = pyo.Param(initialize=config.multiplier)

    unit.outlet_pressure = pyo.Param(initialize=config.outlet_pressure)
    unit.outlet_temperature = pyo.Param(initialize=config.outlet_temperature)

#adding variables and constraints
def add_equations(unit,name,config):
    #local variables
    inlet=unit.control_volume.properties_in[0]
    reservoir_state=unit.control_volume.reservoir_state[0]
    outlet=unit.control_volume.properties_out[0]
    epsilon = 1E-4

    #fix reservoir temperature and pressure
    unit.reservoir_temperature_constraint = pyo.Constraint(
        expr=reservoir_state.temperature==unit.reservoir_temperature
    )
    unit.reservoir_pressure_constraint = pyo.Constraint(
        expr=reservoir_state.pressure==unit.reservoir_pressure
    )
    set_scaling_factor(unit.reservoir_temperature_constraint,1e-2)
    set_scaling_factor(unit.reservoir_pressure_constraint,1e-7)

    #define HCPV
    unit.HCPV = pyo.Var(domain=pyo.NonNegativeReals)

    #create slack variables for initialization
    num_slacks=1
    unit.spos = pyo.Var(range(1,num_slacks+1),domain=pyo.NonNegativeReals,initialize=0)
    unit.sneg = pyo.Var(range(1,num_slacks+1),domain=pyo.NonNegativeReals,initialize=0)
    unit.spos.fix(0)
    unit.sneg.fix(0)
    #create feasibility expression and objective function
    unit.feasibility_expression = pyo.Expression(expr=pyo.quicksum(unit.spos[i]+unit.sneg[i] for i in range(1,num_slacks+1)))
    unit.feasibility_objective = pyo.Objective(expr=unit.feasibility_expression)
    unit.feasibility_objective.deactivate()

    #equality constraints

    #connect inlet state and reservoir state mass balance
    add_state_material_balances(unit,balance_type=idaes.core.MaterialBalanceType.componentTotal,state_1=unit.control_volume.properties_in,state_2=unit.control_volume.reservoir_state,name='material_balance')
    set_scaling_factor(unit.material_balance,1e-1)

    #for now, declare an expression on the average state block to handle liquid viscosity
    reservoir_state.visc_d_phase = pyo.Expression(["Liq"],expr=1e-4)

    #use smooth max so that well pattern can be "turned off" and not get negative flow
    unit.darcys_law = pyo.Constraint(
        expr=inlet.flow_vol==
                #unit.kovr/reservoir_state.visc_d_phase["Liq"]*smooth_max(injection_state.pressure-reservoir_state.pressure,0,epsilon)
                    #+ unit.spos[1]-unit.sneg[1]
                unit.injectivity_index*smooth_max(inlet.pressure-unit.reference_pressure,0,epsilon)*unit.multiplier
                    + unit.spos[1]-unit.sneg[1]
    )
    set_scaling_factor(unit.darcys_law,1e3)

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
            expr=1/unit.Boi*unit.dRfdHCPV*reservoir_state.flow_vol*6.29*unit.correction_factor
            )
    
    #derivative of gas breakthrough curve
    unit.dGbdHCPV = pyo.Expression(
        expr=unit.GB_A*unit.GB_B*unit.HCPV**(unit.GB_B-1)
    )
    #gas breakthrough rate
    unit.q_BRKTH = pyo.Expression(
            expr=unit.dGbdHCPV*reservoir_state.flow_vol*6.29
            )
    
    #outlet state PT
    unit.outlet_temperature_constraint = pyo.Constraint(
        expr=outlet.temperature==unit.outlet_temperature
    )
    unit.outlet_pressure_constraint = pyo.Constraint(
        expr=outlet.pressure==unit.outlet_pressure
    )
    set_scaling_factor(unit.outlet_temperature_constraint,1e-2)
    set_scaling_factor(unit.outlet_pressure_constraint,1e-7)

    #hardcode outlet component flows
    unit.co2_out_constraint = pyo.Constraint(
        expr=outlet.flow_mass_comp['co2']==unit.dGbdHCPV*reservoir_state.flow_mass_comp['co2']
    )
    unit.ch4_out_constraint = pyo.Constraint(
        expr=outlet.flow_mass_comp['ch4']==unit.GOR*unit.q_OIL_PROD+unit.dGbdHCPV*reservoir_state.flow_mass_comp['ch4']
    )
    set_scaling_factor(unit.co2_out_constraint,1e0)
    set_scaling_factor(unit.ch4_out_constraint,1e0)

    #inequality constraints
    #maximum injection pressure
    unit.max_pressure_constraint = pyo.Constraint(
            expr=inlet.pressure<=unit.pressure_max
            )
    #ensure injection pressure is greater than reservoir pressure
    unit.min_pressure_constraint = pyo.Constraint(
            expr=inlet.pressure>=unit.reference_pressure
            )
    unit.min_pressure_constraint.deactivate()
    set_scaling_factor(unit.min_pressure_constraint,1e-7)
    set_scaling_factor(unit.max_pressure_constraint,1e-7)

def guess_scales(unit,name,config):
    inlet=unit.control_volume.properties_in[0]
    reservoir_state=unit.control_volume.reservoir_state[0]
    outlet=unit.control_volume.properties_out[0]

    #variable scales
    set_scaling_factor(inlet.pressure,1e-7)
    set_scaling_factor(inlet.temperature,1e-2)

    set_scaling_factor(reservoir_state.temperature,1e-2)
    set_scaling_factor(reservoir_state.pressure,1e-7)

    set_scaling_factor(unit.outlet_pressure,1e-7)
    set_scaling_factor(unit.outlet_temperature,1e-2)

    #misc scales
    set_scaling_factor(unit.spos[1],1e3)
    set_scaling_factor(unit.sneg[1],1e3)

#define wellpad class
@declare_process_block_class("wellpattern")
class wellpatternData(UnitModelBlockData):
    CONFIG = UnitModelBlockData.CONFIG()
    make_wellpattern_config_block(CONFIG)

    def build(self):
        super(wellpatternData,self).build()
        make_control_volume(self,"control_volume",self.config)
        add_params(self,"params",self.config)
        add_equations(self,"constraints",self.config)
        guess_scales(self,'scales',self.config)
        self.add_port(block=self.control_volume.properties_in,name="inlet")
        self.add_port(block=self.control_volume.properties_out,name="outlet")
    
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

    def custom_propagate_state(self):
        self.control_volume.properties_out[0].temperature.value = pyo.value(self.outlet_temperature)
        self.control_volume.properties_out[0].enth_mol.value = pyo.value(self.control_volume.properties_in[0].enth_mol)
        self.outlet.pressure[0].value = pyo.value(self.outlet_pressure)
        for j in self.control_volume.properties_out.component_list:
            if j == 'co2':
                self.outlet.flow_mol_comp[0,j].value = pyo.value(self.control_volume.properties_in[0].flow_mol)
                self.control_volume.properties_out[0].mole_frac_comp[j].value = 1
            else:
                self.control_volume.properties_out[0].flow_mol_comp[j].value = 0
                self.control_volume.properties_out[0].mole_frac_comp[j].value = 1e-14

    def initialize_states(self):
        self.control_volume.properties_in.initialize()
        self.control_volume.reservoir_state.initialize()
        self.control_volume.properties_out.initialize()

    def initialize(self,solver=None,tee=False,display_after=False):
        print(f'initializing {self.name}')
        self.custom_propagate_state()
        self.initialize_states()
        #activate feasibility problem
        self.activate_feasibility_problem()
        #scale model
        scaled_self = pyo.TransformationFactory('core.scale_model').create_using(self)
        if solver==None:
            solver = pyo.SolverFactory('ipopt')
            solver.options['linear_solver']='ma27'
        res = solver.solve(scaled_self,tee=tee)
        #undo scaling
        pyo.TransformationFactory('core.scale_model').propagate_solution(scaled_self,self)
        if display_after: 
            self.display()
        #deactivate feasibility problem
        self.deactivate_feasibility_problem()
        if res.solver.termination_condition == pyo.TerminationCondition.optimal:
            #create autoscaler
            autoScaler=AutoScaler(overwrite=True)
            autoScaler.scale_variables_by_magnitude(self)
            #autoScaler.scale_constraints_by_jacobian_norm(self)
            print(f'{self.name} initialization solve successful')
        else:
            print(f'{self.name} initialization solve failed, propagating state')
            self.custom_propagate_state()
        return res

    def export_df(self):
        data = {
                "reservoir temperature (K)":pyo.value(self.reservoir_temperature),
                "reservoir pressure (bar)": pyo.value(self.reservoir_pressure)/100000,
                "oil volume factor (STB oil/RB oil)":pyo.value(self.Boi),
                "gas oil ratio (kg gas/STB oil)":pyo.value(self.GOR),
                "production curve fit parameter A":pyo.value(self.PC_A),
                "production curve fit parameter B":pyo.value(self.PC_B),
                "gas breakthrough curve fit parameter A":pyo.value(self.GB_A),
                "gas breakthrough curve fit parameter B":pyo.value(self.GB_B),
                "k overall (m3)":pyo.value(self.kovr),
                "injectivity index (m3)":pyo.value(self.injectivity_index),
                "reference pressure (bar)":pyo.value(self.reference_pressure)/100000,
                "max inlet pressure (bar)":pyo.value(self.pressure_max)/100000,
                "multiplier":pyo.value(self.multiplier),
                "HCPV":pyo.value(self.HCPV),
                "total injection rate (kg/s)":pyo.value(self.control_volume.properties_in[0].flow_mass),
                "total injection rate (mol/s)":pyo.value(self.control_volume.properties_in[0].flow_mol),
                "inlet pressure (bar)":pyo.value(self.inlet.pressure[0])/100000,
                "inlet temperature (K)":pyo.value(self.control_volume.properties_in[0].temperature),
                "correction factor":pyo.value(self.correction_factor),
                "production curve slope":pyo.value(self.dRfdHCPV),
                "gas breakthrough curve slope":pyo.value(self.dGbdHCPV),
                "oil production rate (STB/s)":pyo.value(self.q_OIL_PROD),
                "gas breakthrough rate (RB/s)":pyo.value(self.q_BRKTH),
                "total CO2 out (kg/s)":pyo.value(self.control_volume.properties_out[0].flow_mass_comp['co2']),
                "total ch4 out (kg/s)":pyo.value(self.control_volume.properties_out[0].flow_mass_comp['ch4']),
                "total CO2 out (mol/s)":pyo.value(self.control_volume.properties_out[0].flow_mol_comp['co2']),
                "total ch4 out (mol/s)":pyo.value(self.control_volume.properties_out[0].flow_mol_comp['ch4']),
                }

        return pd.DataFrame(data,index=[self.name])