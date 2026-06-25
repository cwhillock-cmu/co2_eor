#gas pipeline mode

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
from idaes.core.util.initialization import propagate_state
from idaes.core.util.config import is_physical_parameter_block
from idaes.core.util.scaling import set_scaling_factor
from idaes.core.scaling.autoscaling import AutoScaler
from idaes.core.util.math import smooth_abs, safe_log

#config options
def make_config_block(config):
    config.declare(
        "property_package",
        ConfigValue(default=useDefault, domain=is_physical_parameter_block)
    )
    config.declare(
        "has_phase_equilibrium",
        ConfigValue(default=False,domain=In([True,False]))
    )
    config.declare(
        "length",
        ConfigValue(default=500,domain=float)
    )
    config.declare(
        "average_pressure_type",
        ConfigValue(default="nonlinear",domain=In(["constant","linear","nonlinear"]))
    )
    config.declare("average_pressure_weight",
            ConfigValue(default=0.5,domain=float)
            )
    config.declare(
        "allowable_stress",
        ConfigValue(default=2000*100000,domain=float)
    )
    config.declare(
        "property_package_args",
        ConfigBlock(implicit=True)
    )

def make_control_volume(unit,name,config):
    control_volume = ControlVolume0DBlock(
        property_package=config.property_package,
        property_package_args=config.property_package_args,
    )
    setattr(unit,"control_volume",control_volume)
    control_volume.add_state_blocks(has_phase_equilibrium=config.has_phase_equilibrium)
    control_volume.properties_avg=config.property_package.build_state_block(unit.flowsheet().time,defined_state=False,has_phase_equilibrium=config.has_phase_equilibrium)

#add parameters from config
def add_params(unit,config):
    unit.length = pyo.Param(initialize=config.length, units=units.m) #m
    unit.R = pyo.Param(initialize=8.314462, units=units.m**3*units.Pa/units.K/units.mol) #m^3*Pa*K^{-1}*mol^{-1}
    unit.average_pressure_weight = pyo.Param(initialize=config.average_pressure_weight)    
    unit.g = pyo.Param(initialize=9.80665, units=units.m/units.s**2) #m/s^2
    unit.allowable_stress = pyo.Param(initialize=config.allowable_stress, units=units.Pa) #Pa

def add_variables(unit,config):
    unit.diameter = pyo.Var(domain=pyo.NonNegativeReals, bounds=(0.003175,5), initialize=1, units=units.m,)
    unit.roughness = pyo.Var(domain=pyo.NonNegativeReals, bounds=(0,0.1), initialize=0.0018, units=units.m)
    unit.max_pressure = pyo.Var(domain=pyo.NonNegativeReals,bounds=(1,2000*100000),initialize=600*100000, units=units.Pa) #Pa

    unit.area = pyo.Expression(expr=3.1415926*unit.diameter**2/4)

def add_equations(unit,config):
    #local variables for convenience
    inlet = unit.control_volume.properties_in[0]
    outlet = unit.control_volume.properties_out[0]
    average = unit.control_volume.properties_avg[0]
    epsilon = 1e-4

    #create a set of slack variables for feasibility problem
    num_slacks=4
    unit.spos = pyo.Var(range(1,num_slacks+1),domain=pyo.NonNegativeReals,initialize=0)
    unit.sneg = pyo.Var(range(1,num_slacks+1),domain=pyo.NonNegativeReals,initialize=0)
    unit.spos.fix(0)
    unit.sneg.fix(0)
    unit.feasibility_expression = pyo.Expression(expr=pyo.quicksum(unit.spos[i]+unit.sneg[i] for i in range(1,num_slacks+1)))
    unit.feasibility_objective = pyo.Objective(expr=unit.feasibility_expression)
    unit.feasibility_objective.deactivate()

    #create velocity variables
    #inlet velocity 
    inlet.velocity = pyo.Var(domain=pyo.NonNegativeReals,initialize=0,bounds=(0,10), units=units.m/units.s)
    outlet.velocity = pyo.Var(domain=pyo.NonNegativeReals,initialize=0,bounds=(0,10), units=units.m/units.s)

    #mass balance 1
    unit.control_volume.add_phase_component_balances(has_phase_equilibrium=config.has_phase_equilibrium)

    #mass balance 2
    unit.add_state_material_balances(balance_type=idaes.core.MaterialBalanceType.componentPhase,state_1=unit.control_volume.properties_in,state_2=unit.control_volume.properties_avg)
    
    #inlet velocity 3
    unit.inlet_velocity = pyo.Constraint(expr=
        inlet.velocity*unit.area*inlet.dens_mass==inlet.flow_mass
        ) 
    
    #outlet velocity 4
    unit.outlet_velocity = pyo.Constraint(expr=
        outlet.velocity*unit.area*outlet.dens_mass==outlet.flow_mass
        )
    
    #average velocity
    average.velocity = pyo.Expression(expr=
        average.flow_mass/average.dens_mass/unit.area
        )
    
    #average temperature 5
    unit.average_temperature = pyo.Constraint(expr=
        average.temperature*2==inlet.temperature+outlet.temperature
        )
        
    #average pressure 6
    unit.constant_average_pressure = pyo.Constraint(expr=
        average.pressure==inlet.pressure + unit.spos[1]-unit.sneg[1]
        )
    unit.constant_average_pressure.deactivate()
    unit.linear_average_pressure = pyo.Constraint(expr=
        average.pressure==unit.average_pressure_weight*inlet.pressure+(1-unit.average_pressure_weight)*outlet.pressure +unit.spos[1]-unit.sneg[1]
        )
    unit.linear_average_pressure.deactivate()
    unit.nonlinear_average_pressure = pyo.Constraint(expr=
        average.pressure==2/3*(inlet.pressure+outlet.pressure-(inlet.pressure*outlet.pressure)/(inlet.pressure+outlet.pressure)) +unit.spos[1]-unit.sneg[1]
        )
    unit.nonlinear_average_pressure.deactivate()
    if config.average_pressure_type == "nonlinear":
        unit.nonlinear_average_pressure.activate()
    elif config.average_pressure_type == "constant":
        unit.constant_average_pressure.activate()
    else:
        unit.linear_average_pressure.activate()

    #isothermal 7
    unit.isothermal = pyo.Constraint(expr=
        inlet.temperature==outlet.temperature + unit.spos[2] - unit.sneg[2])

    #average reynolds number
    average.Re = pyo.Expression(expr=
        #average.dens_mass*average.velocity*unit.unit.diameter/average.visc_d_phase["Liq"]
        #average.dens_mass*average.velocity*unit.unit.diameter/average.visc_d_phase["Liq"]
        smooth_abs(average.dens_mass*average.velocity*unit.diameter/(15E-6),epsilon)
        #average.dens_mass*average.velocity*unit.diameter/average.visc_d_phase["Liq"]+1
        )
    
    #average friction factor
    average.inverse_f = pyo.Expression(expr=
        4*pyo.log10((unit.roughness/3.7/(unit.diameter)+5.74/((average.Re)**0.9)))**(2) +unit.spos[3]-unit.sneg[3]
        )
    
    #compressibility
    average.Z = pyo.Expression(expr=
        average.pressure/unit.R/average.temperature/average.dens_mass)
    
    #specific gravity
    average.G = pyo.Expression(expr=
        average.dens_mass/1.2754) #standard density of air is 1.2754 kg/m3 (0C 100 kPa)
    
    #general flow equation 8
    unit.GFE = pyo.Constraint(expr=
        #average.flow_mass/1.96==
        #    13.3*(298/101325)*unit.diameter**2.5*((inlet.pressure**2-outlet.pressure**2)*average.inverse_f/(average.G*average.temperature*unit.length*average.Z)+epsilon)**0.5
        (average.flow_mass/1.96)**2==
            (13.3*(298/101325))**2*unit.diameter**5*((inlet.pressure**2-outlet.pressure**2)*average.inverse_f/(average.G*average.temperature*unit.length*average.Z))
        ) #might need to modify constants
    
    #auxiliary constraints
    unit.inlet_pressure_max = pyo.Constraint(
        expr=inlet.pressure<=unit.max_pressure
    )
    unit.outlet_pressure_max = pyo.Constraint(
        expr=outlet.pressure<=unit.max_pressure
    )
    unit.average_pressure_max = pyo.Constraint(
        expr=average.pressure<=unit.max_pressure
    )

    #auxiliary expressions
    unit.Pdrop = pyo.Expression(expr=inlet.pressure-outlet.pressure)
    unit.thickness = pyo.Expression(expr=unit.max_pressure*unit.diameter/(2*(unit.allowable_stress-unit.max_pressure)))

    #misc logic
    unit.misc_logic = pyo.ConstraintList()
    unit.misc_logic.add(expr=inlet.flow_mass>=0)
    unit.misc_logic.add(expr=outlet.flow_mass>=0)
    unit.misc_logic.add(expr=average.flow_mass>=0)
    unit.misc_logic.add(expr=inlet.pressure>=outlet.pressure)
    #unit.misc_logic.deactivate

def guess_scales(unit):
    #local variables
    inlet = unit.control_volume.properties_in[0]
    outlet = unit.control_volume.properties_out[0]
    average = unit.control_volume.properties_avg[0]

    #variable and parameter scaling factors
    set_scaling_factor(unit.diameter,1e2)
    set_scaling_factor(unit.length,1e-3)
    set_scaling_factor(unit.roughness,1e2)
    set_scaling_factor(unit.max_pressure,1e-7)
    set_scaling_factor(unit.allowable_stress,1e-8)
    set_scaling_factor(inlet.flow_mass,1e-2)
    set_scaling_factor(outlet.flow_mass,1e-2)
    set_scaling_factor(average.flow_mass,1e-2)
    set_scaling_factor(inlet.dens_mass,1e-2)
    set_scaling_factor(outlet.dens_mass,1e-2)
    set_scaling_factor(average.dens_mass,1e-2)
    set_scaling_factor(inlet.pressure,1e-7)
    set_scaling_factor(outlet.pressure,1e-7)
    set_scaling_factor(average.pressure,1e-7)
    set_scaling_factor(inlet.temperature,1e-2)
    set_scaling_factor(outlet.temperature,1e-2)
    set_scaling_factor(average.temperature,1e-2)
    set_scaling_factor(average.mw,1e2)

    #equality constraint scaling factors
    set_scaling_factor(unit.control_volume.material_balances,1e-2)
    set_scaling_factor(unit.state_material_balances,1e-2)
    set_scaling_factor(unit.average_temperature,1e-2)
    set_scaling_factor(unit.isothermal,1e-2)
    set_scaling_factor(unit.constant_average_pressure,1e-7)
    set_scaling_factor(unit.linear_average_pressure,1e-7)
    set_scaling_factor(unit.nonlinear_average_pressure,1e-7)
    set_scaling_factor(unit.GFE,1e-1)
    
    #auxiliary constraint scaling factors
    set_scaling_factor(unit.inlet_pressure_max,1e-7)
    set_scaling_factor(unit.outlet_pressure_max,1e-7)
    set_scaling_factor(unit.average_pressure_max,1e-7)

#define gas pipe class
@declare_process_block_class("gasPipe")
class gasPipeData(UnitModelBlockData):
    CONFIG = UnitModelBlockData.CONFIG()
    make_config_block(CONFIG)

    def build(self):
        super(gasPipeData,self).build()
        make_control_volume(self,"control_volume",self.config)

        add_params(self,self.config)
        add_variables(self,self.config)
        add_equations(self,self.config)
        guess_scales(self)
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
        print(f'start initialize function {self.name}')
        #activate feasibility problem
        self.activate_feasibility_problem()
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
            self.print_all()
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
            propagate_state(self.inlet,self.outlet)
        return res
    
    def export_df(self,t=0):
        data = {
            "length (m)":pyo.value(self.length),
            "diameter (m)":pyo.value(self.diameter),
            "roughness (m)":pyo.value(self.roughness),
            "max pressure (bar)":pyo.value(self.max_pressure)/100000,
            "allowable stress (bar)":pyo.value(self.allowable_stress)/100000,
            "thickness (m)":pyo.value(self.thickness),
            "total flowrate (kg/s)":pyo.value(self.control_volume.properties_in[t].flow_mass),
            "pressure drop (bar)":pyo.value(self.Pdrop)/100000,
        }
        for j in self.control_volume.properties_in[t].component_list:
            data.update({f'flowrate {j} (kg/s)':pyo.value(self.control_volume.properties_in[t].flow_mass_comp[j])})
            data.update({f'mass frac {j}':pyo.value(self.control_volume.properties_in[t].mass_frac_comp[j])})
            data.update({f'mole frac {j}':pyo.value(self.control_volume.properties_in[t].mole_frac_comp[j])})
        for state,statename in zip([self.control_volume.properties_in[t],self.control_volume.properties_out[t],self.control_volume.properties_avg[t]],['inlet','outlet','average']):
            data.update({f'{statename} temperature (K)':pyo.value(state.temperature)})
            data.update({f'{statename} pressure (bar)':pyo.value(state.pressure)/100000})
            data.update({f'{statename} density (kg/m3)':pyo.value(state.dens_mass)})
            data.update({f'{statename} velocity (m/s)':pyo.value(state.velocity)})
        data.update({"average Re":pyo.value(self.control_volume.properties_avg[t].Re)})
        data.update({"average 1/f":pyo.value(self.control_volume.properties_avg[t].inverse_f)})
        data.update({"average Z":pyo.value(self.control_volume.properties_avg[t].Z)})
        data.update({"average G":pyo.value(self.control_volume.properties_avg[t].G)})

        return pd.DataFrame(data,index=[self.name])

        