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
            "R",
            ConfigValue(default=8.314462,domain=float)
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

#alternatively, ditch control volume and add states directly
def make_states(unit,name,config):
    unit.inlet = config.property_package.build_state_block(defined_state=False)
    unit.outlet = config.property_package.build_state_block(defined_state=False)
    unit.average = config.property_package.build_state_block(defined_state=False)

#adding parameters from config
def add_params(unit,name,config):
    unit.length = pyo.Param(initialize=config.length) #m
    unit.diameter = pyo.Param(initialize=config.diameter) #m
    unit.roughness = pyo.Param(initialize=config.roughness) #m
    unit.area = pyo.Param(initialize=3.1415926*unit.diameter**2/4) #m^2
    unit.R = pyo.Param(initialize=8.314462) #m^3*Pa*K^{-1}*mol^{-1}

#adding variables and constraints
def add_equations(unit,name,config):
    #create local variables for ease
    inlet = unit.control_volume.properties_in[0]
    outlet = unit.control_volume.properties_out[0]
    average = unit.control_volume.properties_avg
    #inlet=unit.inlet
    #outlet=unit.outlet
    #average=unit.average
    #create a set of slack variables for feasibility problem
    num_slacks=12
    unit.spos = pyo.Var(range(1,num_slacks+1),domain=pyo.NonNegativeReals,initialize=0)
    unit.sneg = pyo.Var(range(1,num_slacks+1),domain=pyo.NonNegativeReals,initialize=0)
    unit.spos.fix(0)
    unit.sneg.fix(0)
    unit.feasibility_expression = pyo.Expression(expr=pyo.quicksum(unit.spos[i]+unit.sneg[i] for i in range(1,num_slacks+1)))
    unit.feasibility_objective = pyo.Objective(expr=unit.feasibility_expression)
    unit.feasibility_objective.deactivate()
    #inlet velocity
    inlet.velocity = pyo.Var(domain=pyo.NonNegativeReals)
    unit.inlet_velocity_constraint = pyo.Constraint(expr=
            inlet.velocity==inlet.flow_mass/inlet.dens_mass/unit.area +unit.spos[1]-unit.sneg[1]
            ) #constraint 1
    #outlet velocity
    outlet.velocity = pyo.Var(domain=pyo.NonNegativeReals)
    unit.outlet_velocity_constraint = pyo.Constraint(expr=
            outlet.velocity==outlet.flow_mass/outlet.dens_mass/unit.area +unit.spos[2]-unit.sneg[2]
            ) #constraint 2
    #average pressure
    unit.average_pressure_constraint = pyo.Constraint(expr=
            average.pressure*3==2*(inlet.pressure**3-outlet.pressure**3)/(inlet.pressure**2-outlet.pressure**2) +unit.spos[3]-unit.sneg[3]
            ) #constraint 3
    #average temperature
    unit.average_temperature_constraint = pyo.Constraint(expr=
            average.temperature*2==inlet.temperature+outlet.temperature +unit.spos[4]-unit.sneg[4]
            ) #constraint 4
    #average mass flow
    unit.average_mass_flow_constraint = pyo.Constraint(expr=
            average.flow_mass*2==inlet.flow_mass+outlet.flow_mass +unit.spos[5]-unit.sneg[5]
            ) #constraint 5
    #average velocity
    average.velocity = pyo.Expression(expr=
            average.flow_mass/average.dens_mass/unit.area 
            )
    #average compressibility
    average.Z = pyo.Expression(expr=
            average.pressure*average.mw/average.dens_mass/unit.R/average.temperature
            )
    #average reynolds
    average.Re = pyo.Expression(expr=
            average.dens_mass*average.velocity*unit.diameter/average.visc_d_phase["Vap"]
            ) 
    #average friction factor
    average.f = pyo.Expression(expr=
            0.25*(pyo.log10(unit.roughness/3.7/unit.diameter+5.74/((average.Re+unit.spos[12]-unit.sneg[12])**0.9)))**(-2)
            )
    #isothermal constraint
    unit.isothermal_constraint=pyo.Constraint(expr=
            inlet.temperature==outlet.temperature +unit.spos[6]-unit.sneg[6]
            ) #constraint 6
    #mass balance
    unit.mass_bal_constraint=pyo.Constraint(expr=
            inlet.flow_mass==outlet.flow_mass +unit.spos[7]-unit.sneg[7]
            ) #constraint 7
    #hydraulic equation
    unit.hydraulic_constraint = pyo.Constraint(expr=
            (inlet.pressure**2-outlet.pressure**2)/(2*average.dens_mass**2*average.velocity**2)
            ==average.Z*unit.R/average.mw*average.temperature*(average.f*unit.length/2/unit.diameter+pyo.log(inlet.pressure/outlet.pressure+unit.spos[11]-unit.sneg[11]))+unit.spos[8]-unit.sneg[8]
            ) #constraint 8
    #useful expression
    unit.deltaP = pyo.Expression(expr=
            inlet.pressure-outlet.pressure
            )
    #define some simplified constraints for initialization
    unit.simple_average_pressure_constraint=pyo.Constraint(expr=
            average.pressure*2==inlet.pressure+outlet.pressure +unit.spos[9]-unit.sneg[9]
            ) #constraint 9
    unit.simple_average_pressure_constraint.deactivate()
    unit.simple_f = pyo.Expression(expr=
            #0.3164*average.Re**(-1/4)
            0.018
            )
    unit.simple_hydraulic_constraint=pyo.Constraint(expr=
            (inlet.pressure**2-outlet.pressure**2)/(2*average.dens_mass**2*average.velocity**2)
            ==average.Z*unit.R/average.mw*average.temperature*unit.simple_f*unit.length/2/unit.diameter +unit.spos[10]-unit.sneg[10]
            ) #constraint 10
    unit.simple_hydraulic_constraint.deactivate()
    
    #inequality constraints
    #pressure bound
    unit.pressure_bound_constraint = pyo.Constraint(
            expr=outlet.pressure<=inlet.pressure
            )

def guess_scales(unit,name,config):
    #create local variables for ease
    inlet = unit.control_volume.properties_in[0]
    outlet = unit.control_volume.properties_out[0]
    average = unit.control_volume.properties_avg
    #inlet=unit.inlet
    #outlet=unit.outlet
    #average=unit.average
    set_scaling_factor(unit.diameter,1e1)
    set_scaling_factor(unit.length,1e-3)
    set_scaling_factor(unit.roughness,1e2)
    set_scaling_factor(unit.area,1e1)
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
    set_scaling_factor(average.mw,1e3)
    set_scaling_factor(average.visc_d_phase["Vap"],1e6)
    set_scaling_factor(unit.average_pressure_constraint,1e-7)
    set_scaling_factor(unit.average_temperature_constraint,1e-2)
    set_scaling_factor(unit.average_mass_flow_constraint,1e-2)
    set_scaling_factor(unit.isothermal_constraint,1e-2)
    set_scaling_factor(unit.mass_bal_constraint,1e-2)
    set_scaling_factor(unit.hydraulic_constraint,1e-3)
    set_scaling_factor(unit.simple_average_pressure_constraint,1e-7)
    set_scaling_factor(unit.simple_hydraulic_constraint,1e-3)
    set_scaling_factor(unit.pressure_bound_constraint,1e-7)

#define pipeline class
@declare_process_block_class("pipeline")
class pipelineData(UnitModelBlockData):
    CONFIG = UnitModelBlockData.CONFIG()
    make_pipeline_config_block(CONFIG)

    def build(self):
        super(pipelineData,self).build()
        make_control_volume(self,"control_volume",self.config)
        #make_states(self,'states',self.config)
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

    def initialize(self,solver='ipopt',tee=False):
        print(f'Start initialization Pipeline')
        self.activate_feasibility_problem()
        self.average_pressure_constraint.deactivate()
        self.hydraulic_constraint.deactivate()
        self.simple_average_pressure_constraint.activate()
        self.simple_hydraulic_constraint.activate()
        #create internal solver
        try:
            local_solver=pyo.SolverFactory(solver)
            local_solver.options['linear_solver']='ma97'
            print(f'loaded passed solver')
        except:
            local_solver=pyo.SolverFactory('ipopt')
            print(f'loaded backup solver')
        #create scaled unit
        scaled_self = pyo.TransformationFactory('core.scale_model').create_using(self)
        res = local_solver.solve(scaled_self,tee=tee)
        self.simple_average_pressure_constraint.deactivate()
        self.simple_hydraulic_constraint.deactivate()
        self.average_pressure_constraint.activate()
        self.hydraulic_constraint.activate()
        res=local_solver.solve(scaled_self,tee=tee)
        #undo scaling
        pyo.TransformationFactory('core.scale_model').propagate_solution(scaled_self,self)
        #create autoscaler
        autoScaler=AutoScaler(overwrite=True)
        autoScaler.scale_variables_by_magnitude(self)
        #autoScaler.scale_constraints_by_jacobian_norm(self)
        self.deactivate_feasibility_problem()
        print(f'End initialization Pipeline')

    def print_expressions(self):
        print(f'Average Velocity:{pyo.value(self.control_volume.properties_avg.velocity)}')
        print(f'Average Re:{pyo.value(self.control_volume.properties_avg.Re)}')
        print(f'Average f:{pyo.value(self.control_volume.properties_avg.f)}')
        print(f'Pressure Drop:{pyo.value(self.deltaP)}')



