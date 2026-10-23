#processing facility surrogate unit model

import idaes
import pyomo.environ as pyo
import pandas as pd
from pyomo.environ import units
from pyomo.common.config import ConfigBlock, ConfigValue, In, ListOf
from idaes.core import (
    ControlVolume0DBlock,
    declare_process_block_class,
    UnitModelBlockData,
    useDefault,
    UnitModelBlock
)
from idaes.core.util.config import is_physical_parameter_block
from idaes.core.util.scaling import set_scaling_factor
from idaes.core.scaling.autoscaling import AutoScaler
from idaes.core.util.math import smooth_abs, safe_log
from idaes.core.util.initialization import propagate_state

from idaes.core.surrogate.surrogate_block import SurrogateBlock

from co2_eor.SSLW import SSLWCosting, SSLWCostingData
from co2_eor.utility.eq_builders import add_pressure_thickness_expression

#config options
def make_config_block(config):
    config.declare(
        "inlet_property_package",
        ConfigValue(default=useDefault,domain=is_physical_parameter_block)
    )
    config.declare(
        "inlet_has_phase_equilibrium",
        ConfigValue(default=False,domain=In([True,False]))
    )
    config.declare(
        "outlet_property_package",
        ConfigValue(default=useDefault,domain=is_physical_parameter_block)
    )
    config.declare(
        "outlet_has_phase_equilibrium",
        ConfigValue(default=False,domain=In([True,False]))
    )
    config.declare(
        "recycle_property_package",
        ConfigValue(default=useDefault,domain=is_physical_parameter_block)
    )
    config.declare(
        "recycle_has_phase_equilibrium",
        ConfigValue(default=False,domain=In([True,False]))
    )
    config.declare(
        "surrogate_input_bounds",
        ConfigValue(default=None,domain=list)
    )
    config.declare(
        "unit_model_parameters",
        ConfigValue(default=None,domain=dict)
    )
    config.declare(
        "global_parameters",
        ConfigValue(default=None,domain=dict)
    )
    config.declare(
        "surrogate_object",
        ConfigValue(default=None)
    )
    config.declare(
        "hyperplane_cuts",
        ConfigValue(default=None,domain=list)
    )
    config.declare(
        "property_package_args",
        ConfigBlock(implicit=True)
    )

#make states for "control volume"
#not using actual idaes control volume block
def make_control_volume(unit,name,config):
    unit.control_volume = pyo.Block()
    unit.control_volume.properties_in = config.inlet_property_package.build_state_block(unit.flowsheet().time,defined_state=True,has_phase_equilibrium=config.inlet_has_phase_equilibrium)
    unit.control_volume.properties_out = config.outlet_property_package.build_state_block(unit.flowsheet().time,defined_state=False,has_phase_equilibrium=config.outlet_has_phase_equilibrium)
    unit.control_volume.properties_recycle = config.recycle_property_package.build_state_block(unit.flowsheet().time,defined_state=True,has_phase_equilibrium=config.recycle_has_phase_equilibrium)
    
def add_global_params(unit,config):
    unit.solvent_MW = pyo.Param(initialize=config.global_parameters["solvent_MW"],units=units.kg/units.mol)
    set_scaling_factor(unit.solvent_MW,1e-2)
    unit.solvent_mass_dens = pyo.Param(initialize=config.global_parameters["solvent_mass_dens"],units=units.kg/units.m**3)
    set_scaling_factor(unit.solvent_mass_dens,1e-2)
    unit.recycle_outlet_temperature = pyo.Param(initialize=config.global_parameters["recycle_outlet_temperature"],units=units.K)
    set_scaling_factor(unit.recycle_outlet_temperature,1e-2)

def add_global_vars(unit,config):
    unit.loop_solvent_flow = pyo.Var(domain=pyo.NonNegativeReals,initialize=1,bounds=(0,10000),units=units.mol/units.s)
    set_scaling_factor(unit.loop_solvent_flow,1e-2)
    unit.makeup_solvent_flow = pyo.Var(domain=pyo.NonNegativeReals,initialize=1,bounds=(0,10),units=units.mol/units.s)
    set_scaling_factor(unit.makeup_solvent_flow,1e2)

def add_units(unit,config):
    #local variables
    inlet = unit.control_volume.properties_in[0]
    outlet = unit.control_volume.properties_out[0]
    recycle = unit.control_volume.properties_recycle[0]

    unit.column = UnitModelBlock(dynamic=False,has_holdup=False)
    unit.column.Nstages = pyo.Var(domain=pyo.NonNegativeReals,initialize=5,bounds=(1,20),units=units.dimensionless)
    unit.column.diameter = pyo.Var(domain=pyo.NonNegativeReals,initialize=1,bounds=(0.1,10),units=units.m)
    unit.column.stage_height = pyo.Param(initialize=config.unit_model_parameters["column"]["stage_height"],units=units.m)
    unit.column.length = pyo.Expression(expr=unit.column.Nstages*unit.column.stage_height)
    unit.column.allowable_stress = pyo.Param(initialize=config.unit_model_parameters["column"]["allowable_stress"],units=units.Pa)
    set_scaling_factor(unit.column.allowable_stress,1e-8)
    unit.column.pressure = pyo.Expression(expr=inlet.pressure)
    set_scaling_factor(unit.column.pressure,1e-7)
    add_pressure_thickness_expression(unit.column,unit.column.pressure,unit.column.diameter,unit.column.allowable_stress)

    unit.compressor1 = UnitModelBlock(dynamic=False,has_holdup=False)
    unit.compressor1.config.declare('compressor',ConfigValue(default=True))
    unit.compressor1.config.declare('thermodynamic_assumption',ConfigValue(default=idaes.models.unit_models.pressure_changer.ThermodynamicAssumption.isentropic))
    unit.compressor1.work_mechanical = pyo.Var(unit.flowsheet().time,domain=pyo.NonNegativeReals,initialize=10000,units=units.W)
    set_scaling_factor(unit.compressor1.work_mechanical,1e-4)

    unit.exchanger1 = UnitModelBlock(dynamic=False,has_holdup=False)
    unit.exchanger1.area = pyo.Var(domain=pyo.NonNegativeReals,initialize=1,units=units.m**2)
    set_scaling_factor(unit.exchanger1.area,1e-1)
    unit.exchanger1.pressure = pyo.Expression(expr=inlet.pressure)
    unit.exchanger1.outlet_temperature = pyo.Var(domain=pyo.NonNegativeReals,bounds=(200,300),initialize=260,units=units.K)
    set_scaling_factor(unit.exchanger1.outlet_temperature,1e-2)
    unit.exchanger1.duty = pyo.Var(domain=pyo.NonNegativeReals,initialize=70000,units=units.W)
    set_scaling_factor(unit.exchanger1.duty,1e-4)

    unit.compressor2 = UnitModelBlock(dynamic=False,has_holdup=False)
    unit.compressor2.config.declare('compressor',ConfigValue(default=True))
    unit.compressor2.config.declare('thermodynamic_assumption',ConfigValue(default=idaes.models.unit_models.pressure_changer.ThermodynamicAssumption.isentropic))
    unit.compressor2.work_mechanical = pyo.Var(unit.flowsheet().time,domain=pyo.NonNegativeReals,initialize=1000000,units=units.W)
    set_scaling_factor(unit.compressor2.work_mechanical,1e-6)
    unit.compressor2.duty = pyo.Var(domain=pyo.Reals,bounds=(None,0),initialize=-200000)
    set_scaling_factor(unit.compressor2.duty,1e-6)

    unit.flash3 = UnitModelBlock(dynamic=False,has_holdup=False)
    unit.flash3.flow_vol = pyo.Var(domain=pyo.NonNegativeReals,initialize=0.1,bounds=(0,10),units=units.m**3/units.s)
    unit.flash3.diameter = pyo.Expression(expr=unit.column.diameter)
    unit.flash3.residence_time = pyo.Param(initialize=config.unit_model_parameters["flash3"]["residence_time"],units=units.s)
    unit.flash3.length = pyo.Var(domain=pyo.NonNegativeReals,initialize=1,bounds=(0,100),units=units.m)
    unit.flash3.length_constraint = pyo.Constraint(expr=
        3.14159*unit.flash3.diameter**2*unit.flash3.length==4*unit.flash3.flow_vol*unit.flash3.residence_time
        )
    unit.flash3.allowable_stress = pyo.Param(initialize=config.unit_model_parameters["flash3"]["allowable_stress"],units=units.Pa)
    unit.flash3.pressure = pyo.Var(domain=pyo.NonNegativeReals,initialize=2e5,bounds=(500,2000e5),units=units.Pa)
    add_pressure_thickness_expression(unit.flash3,unit.flash3.pressure,unit.flash3.diameter,unit.flash3.allowable_stress)
    set_scaling_factor(unit.flash3.flow_vol,1e2)
    set_scaling_factor(unit.flash3.length,1)
    set_scaling_factor(unit.flash3.allowable_stress,1e-8)
    set_scaling_factor(unit.flash3.pressure,1e-5)

    unit.flash1 = UnitModelBlock(dynamic=False,has_holdup=False)
    unit.flash1.flow_vol = pyo.Var(domain=pyo.NonNegativeReals,initialize=0.1,bounds=(0,10),units=units.m**3/units.s)
    unit.flash1.diameter = pyo.Expression(expr=unit.column.diameter)
    unit.flash1.residence_time = pyo.Param(initialize=config.unit_model_parameters["flash1"]["residence_time"],units=units.s)
    unit.flash1.length = pyo.Var(domain=pyo.NonNegativeReals,initialize=1,bounds=(0,100),units=units.m)
    unit.flash1.length_constraint = pyo.Constraint(expr=
        3.14159*unit.flash1.diameter**2*unit.flash1.length==4*unit.flash1.flow_vol*unit.flash1.residence_time
        )
    unit.flash1.allowable_stress = pyo.Param(initialize=config.unit_model_parameters["flash1"]["allowable_stress"],units=units.Pa)
    unit.flash1.pressure = pyo.Expression(expr=inlet.pressure-(inlet.pressure-unit.flash3.pressure)*0.33)
    add_pressure_thickness_expression(unit.flash1,unit.flash1.pressure,unit.flash1.diameter,unit.flash1.allowable_stress)
    set_scaling_factor(unit.flash3.flow_vol,1e2)
    set_scaling_factor(unit.flash3.length,1)
    set_scaling_factor(unit.flash3.allowable_stress,1e-8)

    unit.flash2 = UnitModelBlock(dynamic=False,has_holdup=False)
    unit.flash2.flow_vol = pyo.Var(domain=pyo.NonNegativeReals,initialize=0.1,bounds=(0,10),units=units.m**3/units.s)
    unit.flash2.diameter = pyo.Expression(expr=unit.column.diameter)
    unit.flash2.residence_time = pyo.Param(initialize=config.unit_model_parameters["flash2"]["residence_time"],units=units.s)
    unit.flash2.length = pyo.Var(domain=pyo.NonNegativeReals,initialize=1,bounds=(0,100),units=units.m)
    unit.flash2.length_constraint = pyo.Constraint(expr=
        3.14159*unit.flash2.diameter**2*unit.flash2.length==4*unit.flash2.flow_vol*unit.flash2.residence_time
        )
    unit.flash2.allowable_stress = pyo.Param(initialize=config.unit_model_parameters["flash2"]["allowable_stress"],units=units.Pa)
    unit.flash2.pressure = pyo.Expression(expr=inlet.pressure-(inlet.pressure-unit.flash3.pressure)*0.67)
    add_pressure_thickness_expression(unit.flash2,unit.flash2.pressure,unit.flash2.diameter,unit.flash2.allowable_stress)
    set_scaling_factor(unit.flash3.flow_vol,1e2)
    set_scaling_factor(unit.flash3.length,1)
    set_scaling_factor(unit.flash3.allowable_stress,1e-8)

    unit.solvent_tank = UnitModelBlock(dynamic=False,has_holdup=False)
    unit.solvent_tank.storage_time = pyo.Param(initialize=config.unit_model_parameters["solvent_tank"]["storage_time"],units=units.s)
    unit.solvent_tank.volume = pyo.Expression(expr=unit.makeup_solvent_flow*unit.solvent_MW*unit.solvent_tank.storage_time/unit.solvent_mass_dens)
    set_scaling_factor(unit.solvent_tank.storage_time,1e-6)

    unit.pump1 = UnitModelBlock(dynamic=False,has_holdup=False)
    unit.pump1.config.declare('compressor',ConfigValue(default=True))
    unit.pump1.config.declare('thermodynamic_assumption',ConfigValue(default=idaes.models.unit_models.pressure_changer.ThermodynamicAssumption.pump))
    unit.pump1.control_volume = pyo.Block()
    unit.pump1.control_volume.properties_in = pyo.Block(unit.flowsheet().time)    
    unit.pump1.work_mechanical = pyo.Var(unit.flowsheet().time,domain=pyo.NonNegativeReals,initialize=100000,units=units.W)   
    set_scaling_factor(unit.pump1.work_mechanical,1e-5) 
    unit.pump1.control_volume.properties_in[0].flow_vol = pyo.Var(domain=pyo.NonNegativeReals,initialize=0.1,bounds=(0,10),units=units.m**3/units.s)
    set_scaling_factor(unit.pump1.control_volume.properties_in[0].flow_vol,1e2)
    unit.pump1.deltaP = pyo.Expression(unit.flowsheet().time,expr=inlet.pressure-unit.flash3.pressure)
    unit.pump1.control_volume.properties_in[0].dens_mass = pyo.Expression(expr=unit.solvent_mass_dens)

def add_equations(unit,config):
    #local variables
    inlet = unit.control_volume.properties_in[0]
    outlet = unit.control_volume.properties_out[0]
    recycle = unit.control_volume.properties_recycle[0]

    unit.top_pressure_eq = pyo.Constraint(expr=
        outlet.pressure==inlet.pressure)
    set_scaling_factor(unit.top_pressure_eq,1e-6)

    unit.co2_mass_bal = pyo.Constraint(expr=
        inlet.flow_mol_comp['co2']==outlet.flow_mol_comp['co2']+recycle.flow_mol_comp['co2'])
    unit.ch4_mass_bal = pyo.Constraint(expr=
        inlet.flow_mol_comp['ch4']==outlet.flow_mol_comp['ch4']+recycle.flow_mol_comp['ch4'])
    set_scaling_factor(unit.co2_mass_bal,1e-2)
    set_scaling_factor(unit.ch4_mass_bal,1e-2)

    unit.outlet_temperature_eq = pyo.Constraint(expr=
        recycle.temperature==unit.recycle_outlet_temperature)
    set_scaling_factor(unit.outlet_temperature_eq,1e-2)

    unit.recycle_pressure_eq = pyo.Constraint(expr=
                                              recycle.pressure*recycle.flow_mol
                                              ==recycle.flow_mol_comp['co2']*config.recycle_property_package.co2.pressure_crit
                                              +recycle.flow_mol_comp['ch4']*config.recycle_property_package.ch4.pressure_crit)
    set_scaling_factor(unit.recycle_pressure_eq,1e-6)

    outlet.flow_mol_comp['co2'].lb = 0#-config.surrogate_input_bounds[0][0]
    outlet.flow_mol_comp['ch4'].lb = 0#-config.surrogate_input_bounds[0][1]


def build_surrogate(unit,config):
    inlet = unit.control_volume.properties_in[0]
    outlet = unit.control_volume.properties_out[0]
    recycle = unit.control_volume.properties_recycle[0]

    inputs = [inlet.flow_mol_comp['co2'],inlet.flow_mol_comp['ch4'],inlet.pressure,inlet.temperature,unit.loop_solvent_flow,unit.exchanger1.outlet_temperature,unit.flash3.pressure,unit.column.Nstages]
    outputs = [recycle.flow_mol_comp['co2'],recycle.flow_mol_comp['ch4'],#recycle.pressure,#recycle.temperature,
               outlet.temperature,
               unit.column.diameter,unit.exchanger1.area,unit.exchanger1.duty,unit.pump1.control_volume.properties_in[0].flow_vol,
               unit.pump1.work_mechanical,unit.compressor1.work_mechanical,unit.compressor2.work_mechanical,unit.compressor2.duty,
               unit.flash1.flow_vol,unit.flash2.flow_vol,unit.flash3.flow_vol,unit.makeup_solvent_flow]

    unit.surrogate = SurrogateBlock()
    unit.surrogate.build_model(config.surrogate_object,input_vars=inputs,output_vars=outputs,use_surrogate_bounds=False)

    #build inlet bound constraints
    unit.input_min_cons = pyo.ConstraintList()
    for i,var in enumerate(inputs):
        unit.input_min_cons.add(expr=var>=config.surrogate_input_bounds[1][i])
        set_scaling_factor(unit.input_min_cons[i+1],1/config.surrogate_input_bounds[1][i])
    unit.input_max_cons = pyo.ConstraintList()
    for i,var in enumerate(inputs):
        unit.input_max_cons.add(expr=var<=config.surrogate_input_bounds[0][i])
        set_scaling_factor(unit.input_max_cons[i+1],1/config.surrogate_input_bounds[0][i])

    if config.hyperplane_cuts is not None:
        unit.hyperplane_cuts = pyo.ConstraintList()
        for hp_dict in config.hyperplane_cuts:
            lhs = hp_dict["bias"]
            for i,var in enumerate(inputs):
                lhs += hp_dict["weights"][i]*var
            unit.hyperplane_cuts.add(expr=lhs<=0)

def build_costing_model(unit):
    unit.costing = SSLWCosting()

    unit.column.costing = idaes.core.UnitModelCostingBlock(
        flowsheet_costing_block=unit.costing,
        costing_method=SSLWCostingData.cost_vessel,
        costing_method_arguments={
            "vertical":True,
            "material_type":"Carbon_steel",
            "include_platforms_ladders":True,
            "vessel_diameter":unit.column.diameter,
            "vessel_length":unit.column.length,
            "vessel_thickness":unit.column.thickness,
            "number_of_trays":unit.column.Nstages,
            "tray_material":"CarbonSteel"
        }
    )

    unit.compressor1.costing = idaes.core.UnitModelCostingBlock(
        flowsheet_costing_block=unit.costing,
        costing_method=SSLWCostingData.cost_compressor,
        costing_method_arguments={
            "compressor_type":"Centrifugal",
            "drive_type":"ElectricMotor",
            "material_type":"StainlessSteel"
        }
    )

    unit.exchanger1.costing = idaes.core.UnitModelCostingBlock(
        flowsheet_costing_block=unit.costing,
        costing_method=SSLWCostingData.cost_heat_exchanger,
        costing_method_arguments={
            "hx_type":"Utube",
            "material_type":"StainlessSteelStainlessSteel",
            "tube_length":"12ft",
            "costing_pressure":unit.exchanger1.pressure,
        }
    )

    unit.compressor2.costing = idaes.core.UnitModelCostingBlock(
        flowsheet_costing_block=unit.costing,
        costing_method=SSLWCostingData.cost_compressor,
        costing_method_arguments={
            "compressor_type":"Centrifugal",
            "drive_type":"ElectricMotor",
            "material_type":"StainlessSteel"
        }
    )
    unit.compressor2.costing.number_of_units.fix(5)

    unit.flash3.costing = idaes.core.UnitModelCostingBlock(
        flowsheet_costing_block=unit.costing,
        costing_method=SSLWCostingData.cost_vessel,
        costing_method_arguments={
            "vertical":True,
            "material_type":"Carbon_steel",
            "include_platforms_ladders":False,
            "vessel_diameter":unit.flash3.diameter,
            "vessel_length":unit.flash3.length,
            "vessel_thickness":unit.flash3.thickness,
            "number_of_trays":None
        }
    )

    unit.flash1.costing = idaes.core.UnitModelCostingBlock(
        flowsheet_costing_block=unit.costing,
        costing_method=SSLWCostingData.cost_vessel,
        costing_method_arguments={
            "vertical":True,
            "material_type":"Carbon_steel",
            "include_platforms_ladders":False,
            "vessel_diameter":unit.flash1.diameter,
            "vessel_length":unit.flash1.length,
            "vessel_thickness":unit.flash1.thickness,
            "number_of_trays":None
        }
    )

    unit.flash2.costing = idaes.core.UnitModelCostingBlock(
        flowsheet_costing_block=unit.costing,
        costing_method=SSLWCostingData.cost_vessel,
        costing_method_arguments={
            "vertical":True,
            "material_type":"Carbon_steel",
            "include_platforms_ladders":False,
            "vessel_diameter":unit.flash2.diameter,
            "vessel_length":unit.flash2.length,
            "vessel_thickness":unit.flash2.thickness,
            "number_of_trays":None
        }
    )

    #TODO:implement
    unit.solvent_tank.costing = pyo.Block()

    unit.pump1.costing = idaes.core.UnitModelCostingBlock(
        flowsheet_costing_block=unit.costing,
        costing_method=SSLWCostingData.cost_pump,
        costing_method_arguments={
            "pump_type":"Centrifugal",
            "material_type":"StainlessSteel",
            "pump_type_factor":1.4,
            "motor_type":"open"
        }
    )

    bm_factors = {
        "ST_exchanger":3.17,
        "vertical_vessel":4.16,
        "horizontal_vessel":3.05,
        "pump":3.3,
        "compressor":2.15,
    }

    unit.costing.capital_cost = pyo.Expression(expr=
                                                         unit.column.costing.capital_cost * bm_factors['vertical_vessel']
                                                         + unit.exchanger1.costing.capital_cost * bm_factors["ST_exchanger"]
                                                         + unit.pump1.costing.capital_cost * bm_factors["pump"]
                                                         + (unit.compressor1.costing.capital_cost + unit.compressor2.costing.capital_cost) * bm_factors["compressor"]
                                                         + (unit.flash1.costing.capital_cost + unit.flash2.costing.capital_cost + unit.flash3.costing.capital_cost) * bm_factors["horizontal_vessel"]
    )
    set_scaling_factor(unit.costing.capital_cost,1e-7)

    unit.costing.refrigeration_duty = pyo.Expression(expr=unit.exchanger1.duty)
    unit.costing.electrical_power = pyo.Expression(expr=
                                                   unit.pump1.work_mechanical[0]
                                                   +unit.compressor1.work_mechanical[0]
                                                   +unit.compressor2.work_mechanical[0])

#define processing facility class
@declare_process_block_class("rectisolSurrogate")
class processingFacilityData(UnitModelBlockData):
    CONFIG = UnitModelBlockData.CONFIG()
    make_config_block(CONFIG)

    def build(self):
        super(processingFacilityData,self).build()
        make_control_volume(self,"control_volume",self.config)

        add_global_params(self,self.config)
        add_global_vars(self,self.config)
        add_units(self,self.config)
        add_equations(self,self.config)
        build_surrogate(self,self.config)
        self.add_port(block=self.control_volume.properties_in,name="inlet")
        self.add_port(block=self.control_volume.properties_out,name="outlet")
        self.add_port(block=self.control_volume.properties_recycle,name="recycle")
        build_costing_model(self)

    def initialize_states(self):
        self.control_volume.properties_in.initialize()
        self.control_volume.properties_out.initialize()
        self.control_volume.properties_recycle.initialize()

    def custom_propagate_state(self):
        self.control_volume.properties_recycle[0].pressure.value = 74e5
        self.control_volume.properties_recycle[0].temperature.value = pyo.value(self.recycle_outlet_temperature)
        self.control_volume.properties_out[0].pressure.value = pyo.value(self.control_volume.properties_in[0].pressure)
        self.control_volume.properties_out[0].temperature.value = pyo.value(self.exchanger1.outlet_temperature)

        for j in self.control_volume.properties_in.component_list:
            if j != 'co2':
                self.control_volume.properties_out[0].flow_mol_comp[j].value = self.control_volume.properties_in[0].flow_mol_comp[j].value * 0.9
                self.control_volume.properties_recycle[0].flow_mol_comp[j].value = self.control_volume.properties_in[0].flow_mol_comp[j].value * 0.1
            else:
                self.control_volume.properties_out[0].flow_mol_comp[j].value = self.control_volume.properties_in[0].flow_mol_comp[j].value * 0.5
                self.control_volume.properties_recycle[0].flow_mol_comp[j].value = self.control_volume.properties_in[0].flow_mol_comp[j].value * 0.5

        self.initialize_states()

    def deactivate_surrogate_input_bounds(self):
        self.input_min_cons.deactivate()
        self.input_max_cons.deactivate()

    def activate_surrogate_input_bounds(self):
        self.input_min_cons.activate()
        self.input_max_cons.activate()

    def set_surrogate_input_lbs(self):
        pass
        #TODO: MAKE THIS FUNCTION

    def initialize(self,solver=None,tee=False,display_after=False):
        print(f'initializing {self.name}')
        self.custom_propagate_state()

        #scale model
        scaled_self = pyo.TransformationFactory('core.scale_model').create_using(self)
        if solver==None:
            solver = pyo.SolverFactory('ipopt')
            solver.options['linear_solver']='ma27'
            solver.options['acceptable_tol']=1e-3
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
            self.custom_propagate_state()
        return res
    
    def export_df(self,t=0):
        data = {
        }
        for state,statename in zip([self.control_volume.properties_in[t],self.control_volume.properties_out[t],self.control_volume.properties_recycle[t]],['inlet','outlet','recycle']):
            data.update({f'{statename} total flowrate (kg/s)':pyo.value(state.flow_mass)})
            for j in state.component_list:
                data.update({f'{statename} flowrate {j} (mol/s)':pyo.value(state.flow_mol_comp[j])})
                data.update({f'{statename} mole frac {j}':pyo.value(state.mole_frac_comp[j])})
            data.update({f'{statename} temperature (K)':pyo.value(state.temperature)})
            data.update({f'{statename} pressure (bar)':pyo.value(state.pressure)/100000})
            data.update({f'{statename} density (kg/m3)':pyo.value(state.dens_mass)})
        data.update({f'Column Stages ':pyo.value(self.column.Nstages)})
        data.update({f'Loop Solvent Flow (mol/s)':pyo.value(self.loop_solvent_flow)})
        data.update({f'HX1 Outlet Temperature (K)':pyo.value(self.exchanger1.outlet_temperature)})
        data.update({f'Flash3 Pressure (bar)':pyo.value(self.flash3.pressure)/100000})
        data.update({f'Column Diameter (m)':pyo.value(self.column.diameter)})
        data.update({f'HX1 area (m2)':pyo.value(self.exchanger1.area)})
        data.update({f'HX1 Duty (W)':pyo.value(self.exchanger1.duty)})
        data.update({f'Pump 1 Volflow (m3/s)':pyo.value(self.pump1.control_volume.properties_in[0].flow_vol)})
        data.update({f'Pump 1 Work (W)':pyo.value(self.pump1.work_mechanical[0])})
        data.update({f'Compressor 1 Work (W)':pyo.value(self.compressor1.work_mechanical[0])})
        data.update({f'Compressor 2 Work (W)':pyo.value(self.compressor2.work_mechanical[0])})
        data.update({f'Compressor 2 Duty (W)':pyo.value(self.compressor2.duty)})
        data.update({f'Flash 1 Volflow (m3/s)':pyo.value(self.flash1.flow_vol)})
        data.update({f'Flash 2 Volflow (m3/s)':pyo.value(self.flash2.flow_vol)})
        data.update({f'Flash 3 Volflow (m3/s)':pyo.value(self.flash3.flow_vol)})
        data.update({f'Makeup Solvent (mol/s)':pyo.value(self.makeup_solvent_flow)})
        for blk_name, blk in self.component_map(pyo.Block).items():
            if hasattr(blk,"costing") and hasattr(blk.costing,"capital_cost"):
                data.update({f'{blk_name} capital cost':pyo.value(blk.costing.capital_cost)})
        data.update({f'Total electrical power':pyo.value(self.costing.electrical_power)})
        return pd.DataFrame(data,index=[self.name])
