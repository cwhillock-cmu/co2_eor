#cascade column model adapted from https://doi.org/10.1016/j.compchemeng.2010.02.029

import pyomo.environ as pyo
import pandas as pd
from pyomo.environ import units
from pyomo.common.config import ConfigBlock, ConfigValue, In, ListOf, Bool
import idaes.core

from idaes.core import (
    ControlVolume0DBlock,
    declare_process_block_class,
    UnitModelBlockData,
    useDefault,
    MaterialFlowBasis,
)
from idaes.core.util.config import is_physical_parameter_block
from idaes.core.util.scaling import set_scaling_factor
from idaes.core.scaling.autoscaling import AutoScaler

def make_config_block(config):
    config.declare(
            "vap_property_package",
            ConfigValue(default=useDefault, domain=is_physical_parameter_block),
            )
    config.declare(
            "vap_has_phase_equilibrium",
            ConfigValue(default=False, domain=In([True,False]))
            )
    config.declare(
            "vap_property_package_args",
            ConfigBlock(implicit=True)
            )
    config.declare(
            "liq_property_package",
            ConfigValue(default=useDefault, domain=is_physical_parameter_block),
            )
    config.declare(
            "liq_has_phase_equilibrium",
            ConfigValue(default=False, domain=In([True,False]))
            )
    config.declare(
            "liq_property_package_args",
            ConfigBlock(implicit=True)
            )
    config.declare(
        "solutes",
        ConfigValue(default=None,domain=ListOf(str))
    )
    config.declare(
        "henry_coefficients",
        ConfigValue(default=None,domain=dict)
    )
    config.declare(
        "solvents",
        ConfigValue(default=None,domain=ListOf(str))
    )

def add_states(unit,config):
    unit.control_volume = pyo.Block()
    #add states
    #likely that has_phase_equilibrium flag can be false for each state block
    #manual VLE calculations will be added between the necessary states
    unit.control_volume.bot_in = config.vap_property_package.build_state_block(unit.flowsheet().time,defined_state=True,has_phase_equilibrium=config.vap_has_phase_equilibrium)
    unit.control_volume.bot_out = config.liq_property_package.build_state_block(unit.flowsheet().time,defined_state=False,has_phase_equilibrium=config.liq_has_phase_equilibrium)
    unit.control_volume.top_in = config.liq_property_package.build_state_block(unit.flowsheet().time,defined_state=True,has_phase_equilibrium=config.liq_has_phase_equilibrium)
    unit.control_volume.top_out = config.vap_property_package.build_state_block(unit.flowsheet().time,defined_state=False,has_phase_equilibrium=config.vap_has_phase_equilibrium)

def add_params(unit,config):
    #make_component_lists
    comp_list = []
    for j in list(unit.control_volume.bot_in.component_list):
        if j not in comp_list:
            comp_list.append(j)
    for j in list(unit.control_volume.top_in.component_list):
        if j not in comp_list:
            comp_list.append(j)
    comp_list_pruned = []
    for j in comp_list:
        if j in config.solvents:
            continue
        else:
            comp_list_pruned.append(j)
    unit.comp_list = pyo.Set(initialize=comp_list)
    unit.equilibrium_comps = pyo.Set(initialize=comp_list_pruned)
    unit.solvents = pyo.Set(initialize=config.solvents)
    unit.solutes = pyo.Set(initialize=config.solutes)

    unit.henry_coefficients = pyo.Param(unit.solutes,['A','B'], initialize=config.henry_coefficients)

def add_variables(unit,config):
    unit.num_trays = pyo.Var(domain=pyo.NonNegativeReals,initialize=1)

def add_equations(unit,config):
    #local variables
    bottom_in = unit.control_volume.bot_in[0]
    bottom_out = unit.control_volume.bot_out[0]
    top_in = unit.control_volume.top_in[0]
    top_out = unit.control_volume.top_out[0]

    #component mass balance
    @unit.Constraint(unit.comp_list)
    def comp_mass_balance(b,j):
        expr = 0
        if j in unit.control_volume.bot_in.component_list:
            expr += bottom_in.flow_mol_comp[j] 
            expr -= top_out.flow_mol_comp[j]
        if j in unit.control_volume.top_in.component_list:
            expr += top_in.flow_mol_comp[j]
            expr -= bottom_out.flow_mol_comp[j]
        return expr==0

    #energy balance
    @unit.Constraint()
    def energy_balance(b):
        return sum(bottom_in.get_enthalpy_flow_terms(p) for p in bottom_in.phase_list) + sum(top_in.get_enthalpy_flow_terms(p) for p in top_in.phase_list) == sum(bottom_out.get_enthalpy_flow_terms(p) for p in bottom_out.phase_list) + sum(top_out.get_enthalpy_flow_terms(p) for p in top_out.phase_list)

    ##henrys constants - top and bottom of column
    #unit.H_top = pyo.Expression(unit.solutes)
    #unit.H_bot = pyo.Expression(unit.solutes)

    #henrys law equation
    @unit.Expression(unit.solutes)
    def H_top(b,j):
        return pyo.exp(unit.henry_coefficients[j,'A'] + unit.henry_coefficients[j,'B'] / top_out.temperature)
    @unit.Expression(unit.solutes)
    def H_bot(b,j):
        return pyo.exp(unit.henry_coefficients[j,'A'] + unit.henry_coefficients[j,'B']  / bottom_out.temperature)

    #equilibrium ratios - top and bottom of column
    unit.K_top = pyo.Var(unit.equilibrium_comps,domain=pyo.NonNegativeReals,initialize=0.1)
    unit.K_bot = pyo.Var(unit.equilibrium_comps,domain=pyo.NonNegativeReals,initialize=0.1)

    #equilibrium ratio constraints
    #for solutes use henry's law, for other components use ratio of fugacity
    @unit.Constraint(unit.equilibrium_comps)
    def top_equilibrium_ratio_eq(b,j):
        if j in config.solutes:
            return unit.K_top[j] == unit.H_top[j] / bottom_out.pressure
        else:
            return unit.K_top[j] == top_out.fug_phase_comp('liq',j) / top_out.fug_phase_comp('vap',j) #hard code phase names - don't know how to get liquid over vapor otherwise
    @unit.Constraint(unit.equilibrium_comps)
    def bot_equilibrium_ratio_eq(b,j):
        if j in config.solutes:
            return unit.K_bot[j] == unit.H_bot[j] / bottom_out.pressure
        else:
            return unit.K_bot[j] == bottom_out.fug_phase_comp('liq',j) / bottom_out.fug_phase_comp('vap',j) #hard code phase names - don't know how to get liquid over vapor otherwise

    #ancillary variables
    unit.L_top = pyo.Var(domain=pyo.NonNegativeReals, units=units.mole/units.s) #liquid flow from first tray down
    unit.V_bot = pyo.Var(domain=pyo.NonNegativeReals, units=units.mole/units.s) #vapor flow from first tray up

    #performance expressions
    @unit.Expression(unit.equilibrium_comps)
    def absorption_factor_top(b,j):
        return unit.L_top / (unit.K_top[j]*top_out.flow_mol)
    @unit.Expression(unit.equilibrium_comps)
    def absorption_factor_bot(b,j):
        return bottom_out.flow_mol / (unit.K_bot[j]*unit.V_bot)
    @unit.Expression(unit.equilibrium_comps)
    def stripping_factor_top(b,j):
        return 1 / unit.absorption_factor_top[j]
    @unit.Expression(unit.equilibrium_comps)
    def stripping_factor_bot(b,j):
        return 1 / unit.absorption_factor_bot[j]
    @unit.Expression(unit.equilibrium_comps)
    def absorption_factor_effective(b,j):
        return pyo.sqrt(unit.absorption_factor_bot[j]*(unit.absorption_factor_top[j]+1) + 0.25) - 0.5
    @unit.Expression(unit.equilibrium_comps)
    def stripping_factor_effective(b,j):
        return pyo.sqrt(unit.stripping_factor_top[j]*(unit.stripping_factor_bot[j]+1) + 0.25) - 0.5
    @unit.Expression(unit.equilibrium_comps)
    def recovery_factor_absorption(b,j):
        return (unit.absorption_factor_effective[j] -  1) / (unit.absorption_factor_effective[j]**(unit.num_trays+1) - 1)
    @unit.Expression(unit.equilibrium_comps)
    def recovery_factor_stripping(b,j):
        return (unit.stripping_factor_effective[j] -  1) / (unit.stripping_factor_effective[j]**(unit.num_trays+1) - 1)

    #performance equation
    @unit.Constraint(unit.equilibrium_comps)
    def performance_equation(b,j):
        return top_out.flow_mol_comp[j] == bottom_in.flow_mol_comp[j]*unit.recovery_factor_absorption[j] + top_in.flow_mol_comp[j]*unit.recovery_factor_stripping[j]

    #proposed constraints
    @unit.Constraint()
    def approximate_mass_balance(b):
        return unit.L_top - bottom_out.flow_mol == top_out.flow_mol - unit.V_bot
    @unit.Constraint()
    def vapor_outlet_dew_point(b):
        return sum(top_out.mole_frac_comp[j] / unit.K_top[j] for j in unit.equilibrium_comps) == 1
    @unit.Constraint()
    def liquid_outlet_bubble_point(b):
        return sum(bottom_out.mole_frac_comp[j] * unit.K_bot[j] for j in unit.equilibrium_comps) == 1

    #misc equations
    #isobaric
    @unit.Constraint()
    def top_pressure_bal(b):
        return top_out.pressure == top_in.pressure
    @unit.Constraint()
    def bot_pressure_bal(b):
        return bottom_out.pressure == bottom_in.pressure

    @unit.Constraint()
    def top_temperature_isothermal(b):
        return top_out.temperature == top_in.temperature
    unit.top_temperature_isothermal.deactivate()
    @unit.Constraint()
    def bottom_temperature_isothermal(b):
        return bottom_out.temperature == bottom_in.temperature
    unit.bottom_temperature_isothermal.deactivate()

@declare_process_block_class("cascade")
class cascadeData(UnitModelBlockData):
    CONFIG = UnitModelBlockData.CONFIG()
    make_config_block(CONFIG)

    def build(self):
        super(cascadeData,self).build()
        add_states(self,self.config)
        add_params(self,self.config)
        add_variables(self,self.config)
        add_equations(self,self.config)
        self.add_port(block=self.control_volume.bot_in,name="bottom_inlet")
        self.add_port(block=self.control_volume.bot_out,name="bottom_outlet")
        self.add_port(block=self.control_volume.top_in,name="top_inlet")
        self.add_port(block=self.control_volume.top_out,name="top_outlet")


