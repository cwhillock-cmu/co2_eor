#Pipeline V2 Model

import pyomo.environ as pyo
from pyomo.common.config import ConfigBlock, ConfigValue, In
from idaes.core import (
    ControlVolume0DBlock,
    declare_process_block_class,
    EnergyBalanceType,
    MomentumBalanceType,
    MaterialBalanceType,
    UnitModelBlockData,
    useDefault,
    FlowsheetBlock,
)
from idaes.core.util.config import is_physical_parameter_block
from idaes.core.util.scaling import set_scaling_factor
from idaes.core.scaling.autoscaling import AutoScaler
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.core.scaling.custom_scaler_base import CustomScalerBase
from idaes.core.util.math import safe_log,smooth_max
def safe_log10(a,eps=1e-7):
    return pyo.log10(smooth_max(a,eps,eps=eps))

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
            "diameter",
            ConfigValue(default=0.5,domain=float)
            )
    config.declare(
            "roughness",
            ConfigValue(default=0.0475e-3,domain=float)
            )
    config.declare(
            "alpha", #ambient heat transfer coefficient
            ConfigValue(default=5,domain=float)
            )
    config.declare(
            "ambient_temperature",
            ConfigValue(default=300,domain=float)
            )
    config.declare(
            "average_JT_coefficient",
            ConfigValue(default=0,domain=float)
            )
    config.declare(
            "average_pressure_type",
            ConfigValue(default="constant",domain=In(["constant","linear","nonlinear"]))
            )
    config.declare(
            "heat_balance_type",
            ConfigValue(default="isothermal",domain=In(["isothermal","nonisothermal"]))
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
    unit.length = pyo.Param(initialize=config.length) #m
    unit.diameter = pyo.Param(initialize=config.diameter) #m
    unit.roughness = pyo.Param(initialize=config.roughness) #m
    unit.area = pyo.Param(initialize=3.1415926*unit.diameter**2/4) #m^2
    unit.R = pyo.Param(initialize=8.314462) #m^3*Pa*K^{-1}*mol^{-1}
    unit.alpha = pyo.Param(initialize=config.alpha) #W/m^2K
    unit.ambient_temperature = pyo.Param(initialize=config.ambient_temperature) #K
    unit.JT = pyo.Param(initialize=config.average_JT_coefficient) #K/Pa

#adding variables and constraints
def add_equations(unit,config):
    #create local variables for
    inlet = unit.control_volume.properties_in[0]
    outlet = unit.control_volume.properties_out[0]
    average = unit.control_volume.properties_avg

    #create a set of slack variables for feasibility problem
    num_slacks=5
    unit.spos = pyo.Var(range(1,num_slacks+1),domain=pyo.NonNegativeReals,initialize=0)
    unit.sneg = pyo.Var(range(1,num_slacks+1),domain=pyo.NonNegativeReals,initialize=0)
    unit.spos.fix(0)
    unit.sneg.fix(0)
    unit.feasibility_expression = pyo.Expression(expr=pyo.quicksum(unit.spos[i]+unit.sneg[i] for i in range(1,num_slacks+1)))
    unit.feasibility_objective = pyo.Objective(expr=unit.feasibility_expression)
    unit.feasibility_objective.deactivate()
    
    #create velocity variables
    #inlet velocity 
    inlet.velocity = pyo.Var(domain=pyo.NonNegativeReals)
    outlet.velocity = pyo.Var(domain=pyo.NonNegativeReals)

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
        average.temperature*2==inlet.temperature+outlet.temperature
        )

    #average pressure 6
    unit.constant_average_pressure = pyo.Constraint(expr=
        average.pressure==inlet.pressure
        )
    unit.constant_average_pressure.deactivate()
    unit.linear_average_pressure = pyo.Constraint(expr=
        average.pressure*2==inlet.pressure+outlet.pressure
        )
    unit.linear_average_pressure.deactivate()
    unit.nonlinear_average_pressure = pyo.Constraint(expr=
        average.pressure*3==2*(inlet.pressure**3-outlet.pressure**3)/(inlet.pressure**2-outlet.pressure**2)
        )
    unit.nonlinear_average_pressure.deactivate()
    if config.average_pressure_type == "nonlinear":
        unit.nonlinear_average_pressure.activate()
    elif config.average_pressure_type == "constant":
        unit.constant_average_pressure.activate()
    else:
        unit.linear_average_pressure.activate()

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
        (inlet.PD_ratio+outlet.PD_ratio)/2
        )

    #ambient heat loss
    average.q = pyo.Expression(expr=
        unit.alpha*(unit.ambient_temperature-average.temperature)
        )

    #"difficult" constraints/expressions

    #average reynolds number
    average.Re = pyo.Expression(expr=
        average.dens_mass*average.velocity*unit.diameter/average.visc_d_phase["Vap"]
        )

    #average friction factor
    average.f = pyo.Expression(expr=
        0.25*safe_log10((unit.roughness/3.7/unit.diameter+5.74/((average.Re)**0.9)))**(-2) +unit.spos[2]-unit.sneg[2]
        )

    #hydraulic equation 7
    unit.hydraulic = pyo.Constraint(expr=
        (inlet.pressure**2-outlet.pressure**2)/(2*average.dens_mass**2*average.velocity**2)==average.PD_ratio*average.f/2*unit.length/unit.diameter+average.PD_ratio*safe_log(inlet.pressure/outlet.pressure)-(inlet.PD_ratio-outlet.PD_ratio) +unit.spos[3]-unit.sneg[3]
        )
    
    #heat balance 8
    unit.isothermal = pyo.Constraint(expr=
        inlet.temperature==outlet.temperature+unit.spos[4]-unit.sneg[4]
        )
    unit.isothermal.deactivate()
    unit.nonisothermal = pyo.Constraint(expr=
        inlet.temperature - outlet.temperature == unit.JT*(inlet.pressure-outlet.pressure)-4*average.q/(average.dens_mass*average.velocity*average.cp_mass)*unit.length/unit.diameter +unit.spos[5]-unit.sneg[5]
        )
    unit.nonisothermal.deactivate()
    if config.heat_balance_type == "nonisothermal":
        unit.nonisothermal.activate()
    else:
        unit.isothermal.activate()

    #inequality constraints

    #pressure bound
    unit.outlet_pressure_max = pyo.Constraint(
            expr=outlet.pressure<=inlet.pressure
            )
    unit.outlet_pressure_min = pyo.Constraint(
            expr=outlet.pressure>=74*100000
            )

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
    set_scaling_factor(unit.area,1e1)
    set_scaling_factor(unit.ambient_temperature,1e-2)
    set_scaling_factor(unit.JT,1e7)
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
    set_scaling_factor(average.visc_d_phase["Vap"],1e6)
    set_scaling_factor(average.cp_mass,1e-4)

    #equality constraint scaling factors
    set_scaling_factor(unit.mass_bal,1e-2)
    set_scaling_factor(unit.average_flow,1e-2)
    set_scaling_factor(unit.average_temperature,1e-7)
    set_scaling_factor(unit.constant_average_pressure,1e-7)
    set_scaling_factor(unit.linear_average_pressure,1e-7)
    set_scaling_factor(unit.nonlinear_average_pressure,1e-7)
    set_scaling_factor(unit.hydraulic,1e-3)
    set_scaling_factor(unit.isothermal,1e-2)
    set_scaling_factor(unit.nonisothermal,1e-2)
    
    #inequality constraint scaling factors
    set_scaling_factor(unit.outlet_pressure_max,1e-7)
    set_scaling_factor(unit.outlet_pressure_min,1e-7)

    
#define pipeline class
@declare_process_block_class("pipeline")
class pipelineData(UnitModelBlockData):
    CONFIG = UnitModelBlockData.CONFIG()
    make_pipeline_config_block(CONFIG)

    def build(self):
        super(pipelineData,self).build()
        make_control_volume(self,"control_volume",self.config)

        add_params(self,self.config)
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

    def initialize(self,solver,tee=False,display_after=False):
        print('initializing pipeline')
        #initialize outlet variables as inlet variables
        self.outlet.pressure[0].value = self.inlet.pressure[0].value*0.95
        self.control_volume.properties_avg.pressure.value = self.inlet.pressure[0].value
        self.outlet.temperature[0].value = self.inlet.temperature[0].value
        self.control_volume.properties_avg.temperature.value = self.inlet.temperature[0].value
        self.control_volume.properties_out[0].velocity.value=0
        #activate feasibility problem
        self.activate_feasibility_problem()
        #scale model
        scaled_self = pyo.TransformationFactory('core.scale_model').create_using(self)
        res = solver.solve(scaled_self,tee=tee)
        #undo scaling
        pyo.TransformationFactory('core.scale_model').propagate_solution(scaled_self,self)
        if display_after: 
            self.display()
            self.print_expressions()
        #deactivate feasibility problem
        self.deactivate_feasibility_problem()
        #create autoscaler
        autoScaler=AutoScaler(overwrite=True)
        autoScaler.scale_variables_by_magnitude(self)
        #autoScaler.scale_constraints_by_jacobian_norm(self)
        print('done initializing pipeline')

    def print_expressions(self):
        print(f'Average Velocity:{pyo.value(self.control_volume.properties_avg.velocity)}')
        print(f'Inlet PD ratio:{pyo.value(self.control_volume.properties_in[0].PD_ratio)}')
        print(f'Outlet PD ratio:{pyo.value(self.control_volume.properties_out[0].PD_ratio)}')
        print(f'Average PD ratio:{pyo.value(self.control_volume.properties_avg.PD_ratio)}')
        print(f'Average q:{pyo.value(self.control_volume.properties_avg.q)}')
        print(f'Average Re:{pyo.value(self.control_volume.properties_avg.Re)}')
        print(f'Average f:{pyo.value(self.control_volume.properties_avg.f)}')
        print(f'Pressure Drop:{pyo.value(self.Pdrop)}')

