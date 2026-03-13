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
            ConfigValue(default=300*100000,domain=float)#maximum injection pressure
            )
    config.declare(
            "multiplier",
            ConfigValue(default=1,domain=int) #flow multiplier
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

#alternatively, ditch control volume and add states directly
def make_states(unit,name,config):
    unit.inlet = config.property_package.build_state_block(defined_state=False)
    unit.injection_state = config.property_package.build_state_block(defined_state=False)
    unit.reservoir_state = config.property_package.build_state_block(defined_state=False)

#adding parameters from config
def add_params(unit,name,config):
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
    reservoir_state=unit.control_volume.properties_out[0]

    #define HCPV
    unit.HCPV = pyo.Var(domain=pyo.NonNegativeReals)

    #create slack variables for initialization
    num_slacks=6
    unit.spos = pyo.Var(range(1,num_slacks+1),domain=pyo.NonNegativeReals,initialize=0)
    unit.sneg = pyo.Var(range(1,num_slacks+1),domain=pyo.NonNegativeReals,initialize=0)
    unit.spos.fix(0)
    unit.sneg.fix(0)
    #create feasibility expression and objective function
    unit.feasibility_expression = pyo.Expression(expr=pyo.quicksum(unit.spos[i]+unit.sneg[i] for i in range(1,num_slacks+1)))
    unit.feasibility_objective = pyo.Objective(expr=unit.feasibility_expression)
    unit.feasibility_objective.deactivate()

    #equality constraints
    #mass balance between injection state and inlet state
    unit.mass_bal_constraint1 = pyo.Constraint(
            expr=inlet.flow_mass==injection_state.flow_mass*unit.multiplier +unit.spos[1]-unit.sneg[1]
            ) #constraint 1
    #mass balance between injection state and reservoir state
    unit.mass_bal_constraint2 = pyo.Constraint(
            expr=injection_state.flow_mass==reservoir_state.flow_mass +unit.spos[2]-unit.sneg[2]
            ) #constraint 2
    #isothermal wellbore constraint
    unit.isothermal_constraint = pyo.Constraint(
            expr=inlet.temperature==injection_state.temperature +unit.spos[3]-unit.sneg[3]
            ) #constraint 3
    #for now, BHP is equal to inlet pressure
    unit.pressure_drop_constraint = pyo.Constraint(
            expr=inlet.pressure==injection_state.pressure +unit.spos[4]-unit.sneg[4]
            ) #constraint 4

    #CO2 Injection Rate RB/s
    unit.q_CO2_INJ = pyo.Var(domain=pyo.NonNegativeReals)
    #calculate injection rate
    unit.injection_rate_from_density_constraint = pyo.Constraint(
            expr=unit.q_CO2_INJ==injection_state.flow_mass/injection_state.dens_mass*6.29 +unit.spos[5]-unit.sneg[5]
            ) #constraint 3
    unit.injection_rate_from_darcys_law_constraint = pyo.Constraint(
            expr=unit.q_CO2_INJ==unit.kovr/reservoir_state.visc_d_phase["Vap"]*(injection_state.pressure-reservoir_state.pressure) +unit.spos[6]-unit.sneg[6]
            ) #constraint 4

    #sensitivity curve correction factor
    if config.use_correction_factor:
        unit.correction_factor = pyo.Expression(
                expr=(1-pyo.exp(-unit.SC_B*unit.q_CO2_INJ))/(1-pyo.exp(-unit.SC_B*unit.IR_base))
                )
    else:
        unit.correction_factor = pyo.Expression(expr=1)
    
    #expressions for production rates
    #Oil production rate STB/s
    unit.q_OIL_PROD = pyo.Expression(
            expr=1/unit.Boi*unit.PC_A*unit.PC_B/(unit.HCPV+unit.PC_B)**2*unit.q_CO2_INJ*unit.correction_factor*unit.multiplier
            )
    #gas production rate
    unit.q_GAS_PROD = pyo.Expression(
            expr=unit.GOR*unit.q_OIL_PROD*unit.Boi*unit.multiplier
            )
    #gas breakthrough rate
    unit.q_CO2_BRKTH = pyo.Expression(
            expr=unit.GB_A*unit.GB_B*unit.HCPV**(unit.GB_B-1)*unit.q_CO2_INJ*unit.multiplier
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
    reservoir_state=unit.control_volume.properties_out[0]
    #inlet=unit.inlet
    #reservoir_state=unit.reservoir_state
    set_scaling_factor(unit.PC_A,10)
    set_scaling_factor(unit.PC_B,10)
    set_scaling_factor(unit.GB_A,10)
    set_scaling_factor(unit.kovr,1e13)
    if config.use_correction_factor:
        set_scaling_factor(unit.SC_A,10)
        set_scaling_factor(unit.SC_B,1e-2)
        set_scaling_factor(unit.IR_base,1e3)
    set_scaling_factor(inlet.pressure,1e-7)
    set_scaling_factor(reservoir_state.pressure,1e-7)
    set_scaling_factor(injection_state.pressure,1e-7)
    set_scaling_factor(inlet.temperature,1e-2)
    set_scaling_factor(reservoir_state.temperature,1e-2)
    set_scaling_factor(injection_state.temperature,1e-2)
    set_scaling_factor(inlet.flow_mass,1e1)
    set_scaling_factor(reservoir_state.flow_mass,1e1)
    set_scaling_factor(injection_state.flow_mass,1e1)
    set_scaling_factor(unit.mass_bal_constraint1,1e1)
    set_scaling_factor(unit.mass_bal_constraint2,1e1)
    set_scaling_factor(unit.isothermal_constraint,1e-2)
    set_scaling_factor(unit.pressure_drop_constraint,1e-7)
    set_scaling_factor(unit.q_CO2_INJ,1e3)
    set_scaling_factor(unit.injection_rate_from_density_constraint,1e3)
    set_scaling_factor(unit.injection_rate_from_darcys_law_constraint,1e3)
    set_scaling_factor(unit.BHP_max,1e-7)
    set_scaling_factor(unit.max_BHP_constraint,1e-7)

#define wellpad class
@declare_process_block_class("wellpad")
class wellpadData(UnitModelBlockData):
    CONFIG = UnitModelBlockData.CONFIG()
    make_wellpad_config_block(CONFIG)

    def build(self):
        super(wellpadData,self).build()
        make_control_volume(self,"control_volume",self.config)
        #make_states(self,'states',self.config)
        add_params(self,"params",self.config)
        add_equations(self,"constraints",self.config)
        guess_scales(self,'scales',self.config)
        self.add_inlet_port(block=self.control_volume.properties_in,name="inlet")
        self.add_outlet_port(block=self.control_volume.properties_out,name="reservoir_state")
    
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
        print(f'Start Initialization Wellpad')
        self.activate_feasibility_problem()
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
        res=local_solver.solve(scaled_self,tee=tee)
        #undo scaling
        pyo.TransformationFactory('core.scale_model').propagate_solution(scaled_self,self)
        #self.display()
        #create autoscaler
        autoScaler=AutoScaler(overwrite=True)
        autoScaler.scale_variables_by_magnitude(self)
        #autoScaler.scale_constraints_by_jacobian_norm(self)
        self.deactivate_feasibility_problem()
        print(f'End Initialization Wellpad')

    def print_expressions(self):
        print(f'correction factor={pyo.value(self.correction_factor)}')
        print(f'q_OIL_PROD={pyo.value(self.q_OIL_PROD)}')
        print(f'q_GAS_PROD={pyo.value(self.q_GAS_PROD)}')
        print(f'q_CO2_BRKTH={pyo.value(self.q_CO2_BRKTH)}')

