#Pipeline V2 Model

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
from idaes.core.util.math import smooth_abs

#specify configuration options
def make_pipeline_config_block(config):
    config.declare(
            "property_package",
            ConfigValue(default=useDefault, domain=is_physical_parameter_block),
            )
    config.declare(
            "has_phase_equilibrium",
            ConfigValue(default=False, domain=In([False]))
            )
    config.declare(
            "length",
            ConfigValue(default=20000,domain=float)
            )
    config.declare(
            "alpha", #ambient heat transfer coefficient
            ConfigValue(default=0,domain=float)
            )
    config.declare(
            "ambient_temperature",
            ConfigValue(default=300,domain=float)
            )
    config.declare("average_weight",
            ConfigValue(default=0,domain=float)
            )
    config.declare(
            "height_change",
            ConfigValue(default=0,domain=float)
            )
    config.declare(
            "average_pressure_type",
            ConfigValue(default="nonlinear",domain=In(["constant","linear","nonlinear"]))
            )
    config.declare(
            "heat_balance_type",
            ConfigValue(default="inlet",domain=In(["inlet","nonisothermal","ambient","no_JT"]))
            )
    config.declare(
        "max_pressure",
        ConfigValue(default=600*100000,domain=float)
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
    #create average state block - if want to use time indexing this will need more thought
    control_volume.properties_avg=config.property_package.build_state_block(defined_state=False)    

#adding parameters from config
def add_params(unit,config):
    unit.length = pyo.Param(initialize=config.length, units=units.m) #m
    unit.R = pyo.Param(initialize=8.314462, units=units.m**3*units.Pa/units.K/units.mol) #m^3*Pa*K^{-1}*mol^{-1}
    unit.alpha = pyo.Param(initialize=config.alpha, units=units.W/units.m**2/units.K) #W/m^2K
    unit.ambient_temperature = pyo.Param(initialize=config.ambient_temperature, units=units.K) #K
    unit.average_weight = pyo.Param(initialize=config.average_weight)
    unit.g = pyo.Param(initialize=9.80665, units=units.m/units.s**2) #m/s^2
    unit.height_change = pyo.Param(initialize=config.height_change, units=units.m) #m
    unit.max_pressure = pyo.Param(initialize=config.max_pressure, units=units.Pa) #Pa

def add_variables(unit,config):
    #diameter variable
    unit.diameter = pyo.Var(domain=pyo.NonNegativeReals, bounds=(0,5), initialize=1, units=units.m,)
    unit.roughness = pyo.Var(domain=pyo.NonNegativeReals, bounds=(0,0.1), initialize=0.0018, units=units.m)
    unit.area = pyo.Expression(expr=3.1415926*unit.diameter**2/4)
    
#adding variables and constraints
def add_equations(unit,config):
    #create local variables for
    inlet = unit.control_volume.properties_in[0]
    outlet = unit.control_volume.properties_out[0]
    average = unit.control_volume.properties_avg

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
    inlet.velocity = pyo.Var(domain=pyo.NonNegativeReals, units=units.m/units.s)
    outlet.velocity = pyo.Var(domain=pyo.NonNegativeReals, units=units.m/units.s)

    #"easy" equality constraints/expressions
    
    #mass balance 1
    unit.mass_bal = pyo.Constraint(expr=
        inlet.flow_mass==outlet.flow_mass
        )
    
    #average flow 2
    unit.average_flow = pyo.Constraint(expr=
        average.flow_mass*2==inlet.flow_mass+outlet.flow_mass
        )

    #inlet velocity 3
    unit.inlet_velocity = pyo.Constraint(expr=
        inlet.velocity==inlet.flow_mass/inlet.dens_mass/unit.area
        ) 

    #outlet velocity 4
    unit.outlet_velocity = pyo.Constraint(expr=
        outlet.velocity==outlet.flow_mass/outlet.dens_mass/unit.area
        )

    #average velocity
    average.velocity = pyo.Expression(expr=
        average.flow_mass/average.dens_mass/unit.area
        )

    #average temperature 5
    unit.average_temperature = pyo.Constraint(expr=
        average.temperature==unit.average_weight*inlet.temperature+(1-unit.average_weight)*outlet.temperature
        )

    #inlet pressure density ratio
    inlet.PD_ratio = pyo.Expression(expr=
        inlet.pressure/inlet.dens_mass
        )

    #outlet pressure density ratio
    outlet.PD_ratio = pyo.Expression(expr=
        outlet.pressure/outlet.dens_mass
        )

    #average pressure density ratio
    average.PD_ratio = pyo.Expression(expr=
        average.pressure/average.dens_mass
        #(inlet.PD_ratio+outlet.PD_ratio)/2
        )

    #ambient heat loss
    average.q = pyo.Expression(expr=
        unit.alpha*(unit.ambient_temperature-average.temperature)
        )

    #"difficult" constraints/expressions

    #average pressure 6
    unit.constant_average_pressure = pyo.Constraint(expr=
        average.pressure==inlet.pressure + unit.spos[1]-unit.sneg[1]
        )
    unit.constant_average_pressure.deactivate()
    unit.linear_average_pressure = pyo.Constraint(expr=
        average.pressure==unit.average_weight*inlet.pressure+(1-unit.average_weight)*outlet.pressure +unit.spos[1]-unit.sneg[1]
        )
    unit.linear_average_pressure.deactivate()
    unit.nonlinear_average_pressure = pyo.Constraint(expr=
        average.pressure*3==2*(inlet.pressure**3-outlet.pressure**3)/(inlet.pressure**2-outlet.pressure**2) +unit.spos[1]-unit.sneg[1]
        )
    unit.nonlinear_average_pressure.deactivate()
    if config.average_pressure_type == "nonlinear":
        unit.nonlinear_average_pressure.activate()
    elif config.average_pressure_type == "constant":
        unit.constant_average_pressure.activate()
    else:
        unit.linear_average_pressure.activate()

    #average reynolds number
    average.Re = pyo.Expression(expr=
        smooth_abs(average.dens_mass*average.velocity*unit.diameter/average.visc_d_phase["Liq"])
        )

    #average friction factor
    average.inverse_f = pyo.Expression(expr=
        4*pyo.log10((unit.roughness/3.7/unit.diameter+5.74/((average.Re)**0.9)))**(2) +unit.spos[2]-unit.sneg[2]
        )
    
    #hydraulic equation 7
    unit.hydraulic = pyo.Constraint(expr=
        (inlet.pressure**2-outlet.pressure**2)/2==
            average.pressure*average.dens_mass*average.velocity**2*unit.length/(average.inverse_f*2*unit.diameter)+
                average.pressure*average.dens_mass*average.velocity**2*pyo.log(inlet.pressure/outlet.pressure)-
                    (average.dens_mass*average.velocity)**2*(inlet.PD_ratio-outlet.PD_ratio)+
                        average.pressure*average.dens_mass*unit.g*unit.height_change+
                            unit.spos[3]-unit.sneg[3]
        )
    
    #heat balance 8
    unit.isothermal_inlet = pyo.Constraint(expr=
        inlet.temperature==outlet.temperature+(unit.spos[4]-unit.sneg[4])
        )
    unit.isothermal_inlet.deactivate()
    unit.nonisothermal = pyo.Constraint(expr=
        (outlet.enth_mass - inlet.enth_mass)*average.dens_mass*average.velocity==
            4*average.q*unit.length/unit.diameter +(unit.spos[4]-unit.sneg[4])
        )
    unit.nonisothermal.deactivate()
    unit.isothermal_ambient = pyo.Constraint(expr=
        outlet.temperature==unit.ambient_temperature +(unit.spos[4]-unit.sneg[4])
        )
    unit.isothermal_ambient.deactivate()
    unit.no_JT = pyo.Constraint(expr=
        outlet.temperature==
            unit.ambient_temperature+
                (inlet.temperature-unit.ambient_temperature)*pyo.exp(-4*unit.alpha*unit.length/(average.dens_mass*average.velocity*average.cp_mass*unit.diameter)) +(unit.spos[4]-unit.sneg[4])
        )
    unit.no_JT.deactivate()
    if config.heat_balance_type == "nonisothermal":
        unit.nonisothermal.activate()
    elif config.heat_balance_type =="ambient":
        unit.isothermal_ambient.activate()
    elif config.heat_balance_type =="no_JT":
        unit.no_JT.activate()
    else:
        unit.isothermal_inlet.activate()

    #auxiliary constraints

    #single supercritical phase in pipeline
    unit.outlet_supercritical = pyo.Constraint(
            expr=outlet.temperature_sat>=outlet.temperature_crit
            )
    unit.inlet_supercritical = pyo.Constraint(
            expr=inlet.temperature_sat>=inlet.temperature_crit
            )
    unit.inlet_pressure_max = pyo.Constraint(
        expr=inlet.pressure<=unit.max_pressure
    )
    unit.outlet_pressure_max = pyo.Constraint(
        expr=outlet.pressure<=unit.max_pressure
    )

    #mach number constraint
    average.M = pyo.Expression(expr=average.velocity/average.speed_sound_phase["Liq"])
    unit.mach_number_constraint = pyo.Constraint(
            expr=unit.length<=
                unit.diameter*average.inverse_f/(4*average.heat_capacity_ratio)*
                    ((1-average.M**2)/average.M**2+(1+average.heat_capacity_ratio)/2*
                        pyo.log(average.M**2)/(2/(average.heat_capacity_ratio+1)/(1+(average.heat_capacity_ratio-1)/2*average.M**2)))
            )
    unit.mach_number_constraint.deactivate()

    #auxiliary expressions
    unit.Pdrop = pyo.Expression(expr=inlet.pressure-outlet.pressure)

def guess_scales(unit):
    #create local variables
    inlet = unit.control_volume.properties_in[0]
    outlet = unit.control_volume.properties_out[0]
    average = unit.control_volume.properties_avg
    
    #variable and parameter scaling factors
    set_scaling_factor(unit.diameter,1e1)
    set_scaling_factor(unit.length,1e-3)
    set_scaling_factor(unit.roughness,1e2)
    set_scaling_factor(unit.ambient_temperature,1e-2)
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
    set_scaling_factor(average.visc_d_phase["Liq"],1e6)

    #equality constraint scaling factors
    set_scaling_factor(unit.mass_bal,1e-2)
    set_scaling_factor(unit.average_flow,1e-2)
    set_scaling_factor(unit.average_temperature,1e-2)
    set_scaling_factor(unit.constant_average_pressure,1e-7)
    set_scaling_factor(unit.linear_average_pressure,1e-7)
    set_scaling_factor(unit.nonlinear_average_pressure,1e-7)
    set_scaling_factor(unit.hydraulic,1e-12)
    set_scaling_factor(unit.isothermal_inlet,1e-2)
    set_scaling_factor(unit.nonisothermal,1e-2)
    set_scaling_factor(unit.isothermal_ambient,1e-2)
    set_scaling_factor(unit.no_JT,1e-2)
    
    #auxiliary constraint scaling factors
    set_scaling_factor(unit.inlet_supercritical,1e-2)
    set_scaling_factor(unit.outlet_supercritical,1e-2)
    set_scaling_factor(unit.mach_number_constraint,1e-3)
    set_scaling_factor(unit.inlet_pressure_max,1e-7)
    set_scaling_factor(unit.outlet_pressure_max,1e-7)

    
#define pipeline class
@declare_process_block_class("pipeline")
class pipelineData(UnitModelBlockData):
    CONFIG = UnitModelBlockData.CONFIG()
    make_pipeline_config_block(CONFIG)

    def build(self):
        super(pipelineData,self).build()
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
        print(f'initializing {self.name}')
        #initialize outlet variables as inlet variables
        """
        self.outlet.pressure[0].value = self.inlet.pressure[0].value*0.98
        self.control_volume.properties_avg.pressure.value = self.inlet.pressure[0].value *0.99
        self.outlet.temperature[0].value = self.inlet.temperature[0].value*0.98
        self.control_volume.properties_avg.temperature.value = self.inlet.temperature[0].value *0.99
        """
        self.control_volume.properties_out[0].velocity.value=0
        self.control_volume.properties_avg.velocity.value=0
        #activate feasibility problem
        self.activate_feasibility_problem()
        #scale model
        scaled_self = pyo.TransformationFactory('core.scale_model').create_using(self)
        if solver==None:
            solver = pyo.SolverFactory('ipopt')
            solver.options['linear_solver']='ma97'
            solver.options['tol']=1e-6
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
        print(f'length: {pyo.value(self.length)} m')
        print(f'diameter: {pyo.value(self.diameter)} m')
        print(f'roughness: {pyo.value(self.roughness)} m')
        print(f'heat transfer coefficient: {pyo.value(self.alpha)} W/m2*K')
        print(f'ambient temperature: {pyo.value(self.ambient_temperature)} K')
        print(f'height change: {pyo.value(self.height_change)} m')
        print()

    def print_variables(self):
        print(f'{self.name} variables')
        print(f'Flowrate: {pyo.value(self.inlet.flow_mass[0])} kg/s')
        print(f'Inlet Pressure: {pyo.value(self.inlet.pressure[0])/100000} bar')
        print(f'Inlet Temperature: {pyo.value(self.inlet.temperature[0])} K')
        print(f'Inlet Velocity: {pyo.value(self.control_volume.properties_in[0].velocity)} m/s')
        print(f'Outlet Pressure: {pyo.value(self.outlet.pressure[0])/100000} bar')
        print(f'Outlet Temperature: {pyo.value(self.outlet.temperature[0])} K')
        print(f'Outlet Velocity: {pyo.value(self.control_volume.properties_out[0].velocity)} m/s')
        print(f'Average Pressure: {pyo.value(self.control_volume.properties_avg.pressure)/100000} bar')
        print(f'Average Temperature: {pyo.value(self.control_volume.properties_avg.temperature)} K')
        print(f'Average Viscosity: {pyo.value(self.control_volume.properties_avg.visc_d_phase["Liq"])} Pa*s')
        print()


    def print_expressions(self):
        print(f'{self.name} expressions')
        print(f'Average Velocity: {pyo.value(self.control_volume.properties_avg.velocity)} m/s')
        print(f'Inlet PD ratio: {pyo.value(self.control_volume.properties_in[0].PD_ratio)} J')
        print(f'Outlet PD ratio: {pyo.value(self.control_volume.properties_out[0].PD_ratio)} J')
        print(f'Average PD ratio: {pyo.value(self.control_volume.properties_avg.PD_ratio)} J')
        print(f'Average q: {pyo.value(self.control_volume.properties_avg.q)} W/m2')
        print(f'Average Re: {pyo.value(self.control_volume.properties_avg.Re)}')
        print(f'Average 1/f: {pyo.value(self.control_volume.properties_avg.inverse_f)} 1/J')
        print(f'Average Cp: {pyo.value(self.control_volume.properties_avg.cp_mass)} J/kg')
        print(f'Average Cp/Cv: {pyo.value(self.control_volume.properties_avg.heat_capacity_ratio)}')
        print(f'Average speed of sound: {pyo.value(self.control_volume.properties_avg.speed_sound_phase["Liq"])} m/s')
        print(f'Average M: {pyo.value(self.control_volume.properties_avg.M)}')
        print(f'Pressure Drop: {pyo.value(self.Pdrop)/100000} bar')
        print()

    def print_all(self):
        print()
        self.print_parameters()
        self.print_variables()
        self.print_expressions()

    def export_df(self):
        data = {
                "length (m)":pyo.value(self.length),
                "diameter (m)":pyo.value(self.diameter),
                "roughness (m)":pyo.value(self.roughness),
                "alpha (W/m2/K)":pyo.value(self.alpha),
                "ambient temperature (K)":pyo.value(self.ambient_temperature),
                "height change":pyo.value(self.height_change),
                "flowrate (kg/s)":pyo.value(self.inlet.flow_mass[0]),
                "inlet pressure (bar)":pyo.value(self.inlet.pressure[0])/100000,
                "inlet temperature (K)":pyo.value(self.inlet.temperature[0]),
                "inlet velocity (m/s)":pyo.value(self.control_volume.properties_in[0].velocity),
                "inlet density (kg/m3)":pyo.value(self.control_volume.properties_in[0].dens_mass),
                "outlet pressure (bar)":pyo.value(self.outlet.pressure[0])/100000,
                "outlet temperature (K)":pyo.value(self.outlet.temperature[0]),
                "outlet velocity (m/s)":pyo.value(self.control_volume.properties_out[0].velocity),
                "outlet density (kg/m3)":pyo.value(self.control_volume.properties_out[0].dens_mass),
                "average pressure (bar)":pyo.value(self.control_volume.properties_avg.pressure),
                "average temperature (K)":pyo.value(self.control_volume.properties_avg.temperature),
                "average velocity (m/s)":pyo.value(self.control_volume.properties_avg.velocity),
                "average density (kg/m3)":pyo.value(self.control_volume.properties_avg.dens_mass),
                "average viscosity Pa*s":pyo.value(self.control_volume.properties_avg.visc_d_phase["Liq"]),
                "average q W/m2":pyo.value(self.control_volume.properties_avg.q),
                "average Re":pyo.value(self.control_volume.properties_avg.Re),
                "average 1/f":pyo.value(self.control_volume.properties_avg.inverse_f),
                "average Cp (J/kg)":pyo.value(self.control_volume.properties_avg.cp_mass),
                "average Cp/Cv":pyo.value(self.control_volume.properties_avg.heat_capacity_ratio),
                "average speed of sound (m/s)":pyo.value(self.control_volume.properties_avg.speed_sound_phase["Liq"]),
                "average M":pyo.value(self.control_volume.properties_avg.M),
                "pressure drop (bar)":pyo.value(self.Pdrop)/100000
                }

        return pd.DataFrame(data,index=[self.name])
        


