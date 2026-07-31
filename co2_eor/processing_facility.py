#gas processing facility model

import idaes
import pyomo.environ as pyo
import pandas as pd
from pyomo.environ import units
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
from idaes.core.util.math import smooth_abs, safe_log
from idaes.core.util.initialization import propagate_state

#config options
def make_config_block(config):
    config.declare(
        "property_package",
        ConfigValue(default=useDefault,domain=is_physical_parameter_block)
    )
    config.declare(
        "has_phase_equilibrium",
        ConfigValue(default=False,domain=In([True,False]))
    )
    config.declare(
        "key_component",
        ConfigValue(default='co2',domain=str)
    )
    config.declare(
        "minimum_inlet_pressure",
        ConfigValue(default=5*100000,domain=float)
    )
    config.declare(
        "recycle_pressure",
        ConfigValue(default=2*100000,domain=float)
    )
    config.declare(
        "recycle_temperature",
        ConfigValue(default=298,domain=float)
    )
    config.declare(
        "outlet_pressure",
        ConfigValue(default=2*100000,domain=float)
    )
    config.declare(
        "outlet_temperature",
        ConfigValue(default=298,domain=float)
    )
    config.declare(
        "property_package_args",
        ConfigBlock(implicit=True)
    )

#make states for "control volume"
#not using actual idaes control volume block
def make_control_volume(unit,name,config):
    unit.control_volume = pyo.Block()

    unit.control_volume.properties_in = config.property_package.build_state_block(unit.flowsheet().time,defined_state=True,has_phase_equilibrium=config.has_phase_equilibrium)
    unit.control_volume.properties_out = config.property_package.build_state_block(unit.flowsheet().time,defined_state=False,has_phase_equilibrium=config.has_phase_equilibrium)
    unit.control_volume.properties_recycle = config.property_package.build_state_block(unit.flowsheet().time,defined_state=True,has_phase_equilibrium=config.has_phase_equilibrium)
    
def add_params(unit,config):
    unit.key_component = config.key_component
    unit.minimum_inlet_pressure = config.minimum_inlet_pressure
    unit.recycle_pressure = config.recycle_pressure
    unit.recycle_temperature = config.recycle_temperature
    unit.outlet_pressure = config.outlet_pressure
    unit.outlet_temperature = config.outlet_temperature

def add_equations(unit,config):
    #local variables
    inlet = unit.control_volume.properties_in[0]
    outlet = unit.control_volume.properties_out[0]
    recycle = unit.control_volume.properties_recycle[0]
    comp_list = unit.control_volume.properties_in.component_list

    #mass balance for each component
    @unit.Constraint(comp_list)
    def comp_mass_balance(b,j):
        return inlet.flow_mol_comp[j]==outlet.flow_mol_comp[j]+recycle.flow_mol_comp[j]
    
    #key component recovery equations TBD FOR NOW 
    unit.key_recycle_constraint = pyo.Constraint(expr=
        recycle.mole_frac_comp[unit.key_component]==0.9997
    )
    unit.key_outlet_constraint = pyo.Constraint(expr=
        outlet.mole_frac_comp[unit.key_component]==1-0.9725
        )
    
    #other component recovery equation
    @unit.Constraint(comp_list)
    def other_component_recovery_constraint(b,j):
        if j != unit.key_component:
            return recycle.mole_frac_comp[j]*(1-inlet.mole_frac_comp[unit.key_component])==(1-recycle.mole_frac_comp[unit.key_component])*inlet.mole_frac_comp[j]
        else:
            return pyo.Constraint.Skip
    recycle.mole_frac_comp_eq[unit.key_component].deactivate()

    #auxiliary constraints
    unit.inlet_pressure_constraint = pyo.Constraint(expr=
        inlet.pressure>=unit.minimum_inlet_pressure
    )

    unit.recycle_pressure_constraint = pyo.Constraint(expr=
        recycle.pressure==unit.recycle_pressure
    )

    unit.recycle_temperature_constraint = pyo.Constraint(expr=
        recycle.temperature==unit.recycle_temperature
    )

    unit.outlet_pressure_constraint = pyo.Constraint(expr=
        outlet.pressure==unit.outlet_pressure
    )

    unit.outlet_temperature_constraint = pyo.Constraint(expr=
        outlet.temperature==unit.outlet_temperature
    )

def guess_scales(unit):
    #local variables
    inlet = unit.control_volume.properties_in[0]
    outlet = unit.control_volume.properties_out[0]
    recycle = unit.control_volume.properties_recycle[0]
    comp_list = unit.control_volume.properties_in.component_list

    #parameter scales
    set_scaling_factor(unit.minimum_inlet_pressure,1e-7)
    set_scaling_factor(unit.recycle_pressure,1e-7)
    set_scaling_factor(unit.recycle_temperature,1e-2)
    set_scaling_factor(unit.outlet_pressure,1e-7)
    set_scaling_factor(unit.outlet_temperature,1e-2)

    #constraint scales
    set_scaling_factor(unit.inlet_pressure_constraint,1e-7)
    set_scaling_factor(unit.recycle_pressure_constraint,1e-7)
    set_scaling_factor(unit.recycle_temperature_constraint,1e-2)
    set_scaling_factor(unit.outlet_pressure_constraint,1e-7)
    set_scaling_factor(unit.outlet_temperature_constraint,1e-2)

    set_scaling_factor(unit.key_recycle_constraint,1e1)
    set_scaling_factor(unit.key_outlet_constraint,1e1)
    set_scaling_factor(unit.other_component_recovery_constraint,1e3)

def build_costing_block(unit):
    unit.costing = pyo.Block()
    #simplified processing facility economics
    unit.costing.capex_coefficient = pyo.Param(initialize=42.35) #$/(tonne/yr)
    #unit.costing.capex_coefficient = pyo.Param(initialize=36.5) #$/(tonne/yr)
    unit.costing.oil_price = pyo.Param(initialize=70)

    unit.costing.capital_cost = pyo.Expression(
        expr=unit.costing.capex_coefficient*(unit.control_volume.properties_recycle[0].flow_mass_comp['co2']*3600*24*365/1000)
    )

    unit.costing.opex = pyo.Expression(
        expr=unit.costing.oil_price*0.01*unit.control_volume.properties_in[0].flow_mass_comp['co2']/1.98*35.3147/1000*3600*24*365 #kg/s -> Sm3/s -> Sft3/s -> MSCF/s -> MSCF/yr
    )

#define processing facility class
@declare_process_block_class("processingFacility")
class processingFacilityData(UnitModelBlockData):
    CONFIG = UnitModelBlockData.CONFIG()
    make_config_block(CONFIG)

    def build(self):
        super(processingFacilityData,self).build()
        make_control_volume(self,"control_volume",self.config)

        add_params(self,self.config)
        add_equations(self,self.config)
        guess_scales(self)
        self.add_port(block=self.control_volume.properties_in,name="inlet")
        self.add_port(block=self.control_volume.properties_out,name="outlet")
        self.add_port(block=self.control_volume.properties_recycle,name="recycle")
        build_costing_block(self)

    def initialize(self,solver=None,tee=False,display_after=False):
        print(f'initializing {self.name}')
        self.control_volume.properties_recycle[0].pressure.value = pyo.value(self.recycle_pressure)
        self.control_volume.properties_recycle[0].temperature.value = pyo.value(self.recycle_temperature)
        self.control_volume.properties_out[0].pressure.value = pyo.value(self.outlet_pressure)
        self.control_volume.properties_out[0].temperature.value = pyo.value(self.outlet_temperature)
        for j in self.control_volume.properties_in.component_list:
            if j != self.key_component:
                self.control_volume.properties_out[0].flow_mol_comp[j].value = self.control_volume.properties_in[0].flow_mol_comp[j].value 
                self.control_volume.properties_recycle[0].flow_mol_comp[j].value = 0 
            else:
                self.control_volume.properties_out[0].flow_mol_comp[self.key_component].value = 0
                self.control_volume.properties_recycle[0].flow_mol_comp[self.key_component].value = self.control_volume.properties_in[0].flow_mol_comp[self.key_component].value
        #try this
        #self.control_volume.properties_in.initialize()
        #self.control_volume.properties_in[0].mole_frac_comp.fix()
        #self.control_volume.properties_in[0].mole_frac_comp_eq.deactivate()
        #scale model
        scaled_self = pyo.TransformationFactory('core.scale_model').create_using(self)
        if solver==None:
            solver = pyo.SolverFactory('ipopt')
            solver.options['linear_solver']='ma27'
            solver.options['tol']=1e-6
        res = solver.solve(scaled_self,tee=tee)
        #undo scaling
        pyo.TransformationFactory('core.scale_model').propagate_solution(scaled_self,self)
        if display_after: 
            self.display()
        if res.solver.termination_condition == pyo.TerminationCondition.optimal:
            #create autoscaler
            autoScaler=AutoScaler(overwrite=True)
            autoScaler.scale_variables_by_magnitude(self)
            #autoScaler.scale_constraints_by_jacobian_norm(self)
            print(f'{self.name} initialization solve successful')
        else:
            #self.display()
            #input()
            print(f'{self.name} initialization solve failed, propagating state')
            #manually propagate state in a decent enough way
            for j in self.control_volume.properties_in.component_list:
                if j != self.key_component:
                    self.control_volume.properties_out[0].flow_mol_comp[j].value = self.control_volume.properties_in[0].flow_mol_comp[j].value
                    self.control_volume.properties_recycle[0].flow_mol_comp[j].value = 0
                else:
                    self.control_volume.properties_out[0].flow_mol_comp[self.key_component].value = 0
                    self.control_volume.properties_recycle[0].flow_mol_comp[self.key_component].value = self.control_volume.properties_in[0].flow_mol_comp[self.key_component].value
        return res
    
    def export_df(self,t=0):
        data = {
            "key component":pyo.value(self.key_component),
            "minimum inlet pressure (bar)":pyo.value(self.minimum_inlet_pressure)/100000,
        }
        for state,statename in zip([self.control_volume.properties_in[t],self.control_volume.properties_out[t],self.control_volume.properties_recycle[t]],['inlet','outlet','recycle']):
            data.update({f'{statename} total flowrate (kg/s)':pyo.value(state.flow_mass)})
            for j in state.component_list:
                data.update({f'{statename} flowrate {j} (mol/s)':pyo.value(state.flow_mol_comp[j])})
                data.update({f'{statename} mole frac {j}':pyo.value(state.mole_frac_comp[j])})
            data.update({f'{statename} temperature (K)':pyo.value(state.temperature)})
            data.update({f'{statename} pressure (bar)':pyo.value(state.pressure)/100000})
            data.update({f'{statename} density (kg/m3)':pyo.value(state.dens_mass)})

        return pd.DataFrame(data,index=[self.name])
