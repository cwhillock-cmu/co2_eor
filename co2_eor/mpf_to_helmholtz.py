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
"""
block for converting a modular properties framework EOS stream
into a helmholtz eos stream in IDAES
can convert_all - all mass in MPF stream becomes co2 in helmholtz stream
or convert_co2 - only co2 in MPF stream is co2 in helmholtz stream
"""

#config options
def make_config_block(config):
    config.declare(
        "property_package_in",
        ConfigValue(default=useDefault,domain=is_physical_parameter_block)
    )
    config.declare(
        "property_package_out",
        ConfigValue(default=useDefault,domain=is_physical_parameter_block)
    )
    config.declare(
        "has_phase_equilibrium",
        ConfigValue(default=False,domain=bool)
    )
    config.declare(
        "conversion_type",
        ConfigValue(default='convert_all_mol',domain=In(['convert_all_mass','convert_co2','convert_all_mol']))
    )

def make_control_volume(unit,name,config):
    unit.control_volume = pyo.Block()

    unit.control_volume.properties_in = config.property_package_in.build_state_block(unit.flowsheet().time,defined_state=True,has_phase_equilibrium=config.has_phase_equilibrium)
    unit.control_volume.properties_out = config.property_package_out.build_state_block(unit.flowsheet().time,defined_state=False,has_phase_equilibrium=config.has_phase_equilibrium)
    
def add_equations(unit,config):
    #local variables
    inlet = unit.control_volume.properties_in[0]
    outlet = unit.control_volume.properties_out[0]

    #temperature equality
    unit.temperature_equality = pyo.Constraint(expr=
        inlet.temperature==outlet.temperature
        )
    #pressure equality
    unit.pressure_equality = pyo.Constraint(expr=
        inlet.pressure==outlet.pressure
        )
    
    #convert mass
    if config.conversion_type == 'convert_all_mass':
        unit.mass_constraint = pyo.Constraint(expr=
            inlet.flow_mass == outlet.flow_mass
            )
    elif config.conversion_type == 'convert_co2':
        unit.mass_constraint = pyo.Constraint(expr=
            inlet.flow_mass_comp['co2'] == outlet.flow_mass
            )
    elif config.conversion_type == 'convert_all_mol':
        unit.mass_constraint = pyo.Constraint(expr=
            inlet.flow_mol == outlet.flow_mol
            )
    else:
        raise ValueError("invalid conversion type")
    
def guess_scales(unit):
    #local variables
    inlet = unit.control_volume.properties_in[0]
    outlet = unit.control_volume.properties_out[0]

    #variable scales
    set_scaling_factor(inlet.temperature,1e-2)
    set_scaling_factor(inlet.pressure,1e-7)
    set_scaling_factor(outlet.temperature,1e-2)
    set_scaling_factor(outlet.pressure,1e-7)

    #constraint scales
    set_scaling_factor(unit.temperature_equality,1e-2)
    set_scaling_factor(unit.pressure_equality,1e-2)
    set_scaling_factor(unit.mass_constraint,1e1)

#define converter unit class
@declare_process_block_class("mpf_helmholtz_converter")
class mpf_helmholtz_converterData(UnitModelBlockData):
    CONFIG = UnitModelBlockData.CONFIG()
    make_config_block(CONFIG)

    def build(self):
        super(mpf_helmholtz_converterData,self).build()
        make_control_volume(self,"control_volume",self.config)

        add_equations(self,self.config)
        guess_scales(self)
        self.add_inlet_port(block=self.control_volume.properties_in,name="inlet")
        self.add_outlet_port(block=self.control_volume.properties_out,name="outlet")

    def initialize(self,solver=None,tee=False,display_after=False):
        print(f'initializing {self.name}')
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
            print(f'{self.name} initialization solve failed, propagating state')
            propagate_state(self.inlet,self.outlet)
        print(f'{self.name} initialization complete')
        return res
    
    def export_df(self,t=0):
        data = {
            "conversion_type":self.config.conversion_type,
        }
        for state,statename in zip([self.control_volume.properties_in[t],self.control_volume.properties_out[t]],['inlet','outlet']):
            data.update({f'{statename} total flowrate (kg/s)':pyo.value(state.flow_mass)})
            for j in state.component_list:
                data.update({f'{statename} flowrate {j} (mol/s)':pyo.value(state.flow_mol_comp[j])})
                data.update({f'{statename} mole frac {j}':pyo.value(state.mole_frac_comp[j])})
            data.update({f'{statename} temperature (K)':pyo.value(state.temperature)})
            data.update({f'{statename} pressure (bar)':pyo.value(state.pressure)/100000})

        return pd.DataFrame(data,index=[self.name])
