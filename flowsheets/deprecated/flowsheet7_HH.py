#test flowsheet

import pyomo.environ as pyo
import idaes.core
from pyomo.network import Arc
from pyomo.environ import units
import contextlib
import numpy as np
import sys
import readchar
import idaes.logger as idaeslog

import idaes.models.properties.general_helmholtz as idaesHelmholtz
import idaes.models.unit_models.pressure_changer as idaesPressureChanger

from idaes.core import MomentumBalanceType
from co2_eor.splitter_unit import EnergySplittingType
from co2_eor.mixer_unit import MomentumMixingType
from idaes.core.util.scaling import set_scaling_factor

from co2_eor import mixer,splitter, heater

from idaes.models.properties.modular_properties.base.generic_property import GenericParameterBlock
from co2_eor.liq_pipe import liqPipe
from co2_eor.wellpattern_unit_HH import wellpattern_HH
import co2_eor.MPF.thermo_config as thermo_config
from co2_eor.gas_pipe import gasPipe
from co2_eor.processing_facility import processingFacility
#from co2_eor.surrogates.rectisol_surrogate import rectisolSurrogate
from co2_eor.mpf_to_helmholtz import mpf_helmholtz_converter

from co2_eor.SSLW import SSLWCosting, SSLWCostingData
from co2_eor.SSLW import CompressorType, CompressorDriveType, CompressorMaterial
from co2_eor.SSLW import VesselMaterial

from co2_eor.utility.solvers import ipopt, conopt

m = pyo.ConcreteModel()
m.fs = idaes.core.FlowsheetBlock(dynamic=False)
m.fs.props_helmholtz = idaesHelmholtz.HelmholtzParameterBlock(
        pure_component="CO2",amount_basis=idaesHelmholtz.AmountBasis.MOLE,
        state_vars=idaesHelmholtz.StateVars.PH,phase_presentation=idaesHelmholtz.PhaseType.MIX,
        #has_phase_equilibrium=False,
        )
m.fs.props_mix_vap_cubic = GenericParameterBlock(**thermo_config.configuration_vap_cubic)
m.fs.props_mix_liq_cubic = GenericParameterBlock(**thermo_config.configuration_liq_cubic)

m.fs.costing = SSLWCosting()

#define default unit configurations 
compressor_config = {
        "property_package": m.fs.props_helmholtz,
        "dynamic":False,
        "compressor":True,
        "thermodynamic_assumption":idaesPressureChanger.ThermodynamicAssumption.isentropic,
        }

liqPipe_config = {
        "property_package":m.fs.props_helmholtz,
        "alpha":5,
        "ambient_temperature":293.15,
        "average_pressure_type":'linear',
        "heat_balance_type":'nonisothermal',
        "average_pressure_weight":0.5,
        "average_temperature_weight":0.0,
        "allowable_stress":1380*100000,
        "height_change":0,
        }

mixer_config = {
        "momentum_mixing_type":MomentumMixingType.inequality
        }

splitter_config = {
        "property_package":m.fs.props_helmholtz,
        "ideal_separation":False,
        "energy_split_basis":EnergySplittingType.equal_molar_enthalpy,
        "momentum_balance_type":MomentumBalanceType.none,
        }

wellpattern_config = {
    "primary_property_package":m.fs.props_helmholtz,
    "secondary_property_package":m.fs.props_mix_vap_cubic,
    "temperature":300,
    "pressure":413*100000, #413
    "Boi":1.3,
    "GOR":84.42,
    "PC_A":0.4,
    "PC_B":0.6,
    "GB_A":0.62,
    "GB_B":1.12,
    "kovr":2.93E-14,
    "injectivity_index":2.93e-10,
    "reference_pressure":413*100000,
    "use_correction_factor":False,
    "depth":1000,
    "pressure_max":650*100000, #650?
    "outlet_pressure":30*100000,
    "outlet_temperature":300,
    "multiplier":3,
}

heater_config = {
    "property_package":m.fs.props_helmholtz,
}

gasPipe_config = {
    "property_package":m.fs.props_mix_vap_cubic,
    "average_pressure_type":'nonlinear',
    "allowable_stress":1380*100000,
}

liqpipe_costing_config = {
    "flowsheet_costing_block":m.fs.costing,
    "costing_method":SSLWCostingData.cost_horizontal_vessel,
    "costing_method_arguments":{
        "material_type":VesselMaterial.StainlessSteel316,
        "include_platforms_ladders":False,
    }
}

gaspipe_costing_config = {
    "flowsheet_costing_block":m.fs.costing,
    "costing_method":SSLWCostingData.cost_horizontal_vessel,
    "costing_method_arguments":{
        "material_type":VesselMaterial.StainlessSteel316,
        "include_platforms_ladders":False,
    }
}

compressor_costing_config={
    "flowsheet_costing_block":m.fs.costing,
    "costing_method":SSLWCostingData.cost_compressor,
    "costing_method_arguments":{
        "compressor_type":CompressorType.Centrifugal,
        "drive_type":CompressorDriveType.ElectricMotor,
        "material_type":CompressorMaterial.StainlessSteel,
    }
}
#build units
m.fs.pipe0 = liqPipe(**liqPipe_config,length=15000)
#m.fs.pipe0.diameter.fix(5)
m.fs.pipe0.diameter = 1
m.fs.pipe0.roughness.fix(0.0475e-3)
m.fs.pipe0.max_pressure.fix(140*100000)
m.fs.pipe0.costing = idaes.core.UnitModelCostingBlock(**liqpipe_costing_config)

m.fs.mix0 = mixer(**mixer_config,property_package=m.fs.props_helmholtz,inlet_list=['from_pipe0','from_purge_splitter'])

m.fs.main_comp = idaesPressureChanger.PressureChanger(**compressor_config)
m.fs.main_comp.efficiency_isentropic.fix(0.85)
m.fs.main_comp.costing = idaes.core.UnitModelCostingBlock(**compressor_costing_config)

m.fs.chiller = heater(**heater_config)

m.fs.pipe1 = liqPipe(**liqPipe_config,length=2500)
#m.fs.pipe1.diameter.fix(5)
m.fs.pipe1.diameter = 1
m.fs.pipe1.roughness.fix(0.0475e-3)
m.fs.pipe1.max_pressure.fix(650*100000)
m.fs.pipe1.costing = idaes.core.UnitModelCostingBlock(**liqpipe_costing_config)

m.fs.split1 = splitter(**splitter_config,outlet_list=['to_well1','to_pipe2'])

m.fs.well1 = wellpattern_HH(**wellpattern_config)
m.fs.well1.HCPV.fix(0.8)

m.fs.pipe2 = liqPipe(**liqPipe_config,length=3000)
#m.fs.pipe2.diameter.fix(5)
m.fs.pipe2.diameter = 1
m.fs.pipe2.roughness.fix(0.0475e-3)
m.fs.pipe2.max_pressure.fix(650*100000)
m.fs.pipe2.costing = idaes.core.UnitModelCostingBlock(**liqpipe_costing_config)

m.fs.well2 = wellpattern_HH(**wellpattern_config)
m.fs.well2.HCPV.fix(0.7)

m.fs.pipe3 = gasPipe(**gasPipe_config,length=3000)
#m.fs.pipe3.diameter.fix(5)
m.fs.pipe3.diameter = 1
m.fs.pipe3.roughness.fix(0.0475e-3)
m.fs.pipe3.max_pressure.fix(50*100000)
m.fs.pipe3.costing = idaes.core.UnitModelCostingBlock(**gaspipe_costing_config)

m.fs.mix1 = mixer(**mixer_config,property_package=m.fs.props_mix_vap_cubic,inlet_list=['from_well1','from_pipe3'])

m.fs.pipe4 = gasPipe(**gasPipe_config,length=2600)
#m.fs.pipe4.diameter.fix(5)
m.fs.pipe4.diameter = 1
m.fs.pipe4.roughness.fix(0.0473e-3)
m.fs.pipe4.max_pressure.fix(50*100000)
m.fs.pipe4.costing = idaes.core.UnitModelCostingBlock(**gaspipe_costing_config)

m.fs.processing_facility=processingFacility(
    inlet_property_package=m.fs.props_mix_vap_cubic,
    outlet_property_package=m.fs.props_mix_vap_cubic,
    recycle_property_package=m.fs.props_mix_vap_cubic,
    key_component='co2',
    minimum_inlet_pressure=10*100000,
    recycle_pressure=74*100000,
    recycle_temperature=350,
    outlet_pressure=5*100000
)

m.fs.eos_converter = mpf_helmholtz_converter(
    property_package_in = m.fs.props_mix_liq_cubic,
    property_package_out = m.fs.props_helmholtz,
    conversion_type = 'convert_co2'
)

m.fs.recycle_comp = idaesPressureChanger.PressureChanger(**compressor_config)
m.fs.recycle_comp.efficiency_isentropic.fix(0.85)
m.fs.recycle_comp.costing = idaes.core.UnitModelCostingBlock(**compressor_costing_config)

m.fs.purge_splitter = splitter(**splitter_config,outlet_list=['to_purge','to_mix0'])
#m.fs.purge_splitter.positive_purge_eq1 = pyo.Constraint(expr=m.fs.purge_splitter.to_purge.flow_mol_comp[0,'co2']>=0) #enforce material does not go into system from purge stream
#m.fs.purge_splitter.positive_purge_eq2 = pyo.Constraint(expr=m.fs.purge_splitter.to_purge.flow_mol_comp[0,'ch4']>=0)
m.fs.purge_splitter.to_purge.flow_mol[0].lb =0
m.fs.purge_splitter.to_purge.pressure[0].fix(74*100000)

#create streams
m.fs.s_pipe0_mix0 = Arc(source=m.fs.pipe0.outlet,destination=m.fs.mix0.from_pipe0)
m.fs.s_mix0_main_comp = Arc(source=m.fs.mix0.outlet,destination=m.fs.main_comp.inlet)
m.fs.s_main_comp_chiller = Arc(source=m.fs.main_comp.outlet,destination=m.fs.chiller.inlet)
m.fs.s_chiller_pipe1 = Arc(source=m.fs.chiller.outlet,destination=m.fs.pipe1.inlet)
m.fs.s_pipe1_split1 = Arc(source=m.fs.pipe1.outlet,destination=m.fs.split1.inlet)
m.fs.s_split1_well1 = Arc(source=m.fs.split1.to_well1,destination=m.fs.well1.inlet)
m.fs.s_well1_mix1 = Arc(source=m.fs.well1.outlet,destination=m.fs.mix1.from_well1)
m.fs.s_split1_pipe2 = Arc(source=m.fs.split1.to_pipe2,destination=m.fs.pipe2.inlet)
m.fs.s_pipe2_well2 = Arc(source=m.fs.pipe2.outlet,destination=m.fs.well2.inlet)
m.fs.s_well2_pipe3 = Arc(source=m.fs.well2.outlet,destination=m.fs.pipe3.inlet)
m.fs.s_pipe3_mix1 = Arc(source=m.fs.pipe3.outlet,destination=m.fs.mix1.from_pipe3)
m.fs.s_mix1_pipe4 = Arc(source=m.fs.mix1.outlet,destination=m.fs.pipe4.inlet)
m.fs.s_pipe4_processing_facility = Arc(source=m.fs.pipe4.outlet,destination=m.fs.processing_facility.inlet)
#m.fs.s_pipe4_preheater = Arc(source=m.fs.pipe4.outlet,destination=m.fs.preheater.inlet)
#m.fs.s_preheater_processing_facility = Arc(source=m.fs.preheater.outlet,destination=m.fs.processing_facility.inlet)
m.fs.s_processing_facility_eos_converter = Arc(source=m.fs.processing_facility.recycle,destination=m.fs.eos_converter.inlet)
#m.fs.s_pipe4_precompressor = Arc(source=m.fs.pipe4.outlet,destination=m.fs.precompressor.inlet)
#m.fs.s_precompressor_processing_facility = Arc(source=m.fs.precompressor.outlet,destination=m.fs.processing_facility.inlet)
#m.fs.s_pipe4_processing_facility = Arc(source=m.fs.pipe4.outlet,destination=m.fs.processing_facility.inlet)
#m.fs.s_processing_facility_recycle_comp = Arc(source=m.fs.processing_facility.recycle,destination=m.fs.recycle_comp.inlet)
m.fs.s_eos_converter_recycle_comp = Arc(source=m.fs.eos_converter.outlet,destination=m.fs.recycle_comp.inlet)
m.fs.s_recycle_comp_purge_splitter = Arc(source=m.fs.recycle_comp.outlet,destination=m.fs.purge_splitter.inlet)
m.fs.s_purge_splitter_mix0 = Arc(source=m.fs.purge_splitter.to_mix0,destination=m.fs.mix0.from_purge_splitter)

pyo.TransformationFactory("network.expand_arcs").apply_to(m)

#fix degrees of freedom
m.fs.pipe0.inlet.pressure[0].fix(100*100000)
m.fs.pipe0.inlet.enth_mol[0].fix(m.fs.props_helmholtz.htpx(T=300*units.K,p=100*100000*units.Pa,amount_basis=idaesHelmholtz.AmountBasis.MOLE))
#m.fs.pipe0.inlet.temperature[0].fix(300)

m.fs.main_comp.outlet.pressure[0].fix(650*100000)
m.fs.chiller.outlet_temp_eq = pyo.Constraint(expr=m.fs.chiller.control_volume.properties_out[0].temperature==318)
#m.fs.chiller.control_volume.properties_out[0].temperature.fix(310)

m.fs.recycle_comp.outlet.pressure[0].fix(86*100000)

m.fs.purge_splitter.split_fraction[0,'to_purge'].fix(1e-14)
m.fs.purge_splitter.to_mix0.pressure[0].fix(85*100000)

#m.fs.well1.inlet.pressure[0].fix(450*100000)
#m.fs.well2.inlet.pressure[0].fix(475*100000)

#fix pressure drops instead of diameters
m.fs.pipe0.Pdrop_constraint = pyo.Constraint(expr=m.fs.pipe0.Pdrop==1*100000)
m.fs.pipe1.Pdrop_constraint = pyo.Constraint(expr=m.fs.pipe1.Pdrop==1*100000)
m.fs.pipe2.Pdrop_constraint = pyo.Constraint(expr=m.fs.pipe2.Pdrop==1*100000)
m.fs.pipe3.Pdrop_constraint = pyo.Constraint(expr=m.fs.pipe3.Pdrop==16*100000)
m.fs.pipe4.Pdrop_constraint = pyo.Constraint(expr=m.fs.pipe4.Pdrop==2.5*100000)
set_scaling_factor(m.fs.pipe0.Pdrop_constraint,1e-6)
set_scaling_factor(m.fs.pipe1.Pdrop_constraint,1e-6)
set_scaling_factor(m.fs.pipe2.Pdrop_constraint,1e-6)
set_scaling_factor(m.fs.pipe3.Pdrop_constraint,1e-6)
set_scaling_factor(m.fs.pipe4.Pdrop_constraint,1e-6)

m.fs.split1.to_well1.pressure[0].fix(600*100000)
m.fs.split1.to_pipe2.pressure[0].fix(600*100000)

m.fs.mix0.outlet.pressure[0].fix(84*100000)
m.fs.mix1.outlet.pressure[0].fix(13.5*100000)

#ancillary degrees of freedom?
#m.fs.split1.to_well1.pressure[0].fix(450*100000)
#m.fs.split1.to_pipe2.pressure[0].fix(530*100000)

#initialization degrees of freedom
m.fs.pipe0.inlet.flow_mol[0].fix(160)

print(f'starting solve')

def initialization_type1(m):
    from pyomo.network import SequentialDecomposition
    from idaes.core.util.initialization import propagate_state
    seq = SequentialDecomposition()
    seq.options.select_tear_method = "heuristic"
    seq.options.tear_method = "Wegstein"
    seq.options.iterLim = 15

    G = seq.create_graph(m)
    heauristic_tear_set = seq.tear_set_arcs(G,method="heuristic")
    order = seq.calculation_order(G)

    print("tear set")
    for o in heauristic_tear_set:
        print(o.name)
    print("calculation order")
    for o in order:
        print(f'{o[0].name}')

    #input('paused')

    tear_guesses = {
        "flow_mol":{
            0:600
        },
        #"temperature":{0:300},
        "enth_mol":{0:m.fs.props_helmholtz.htpx(T=300*units.K,p=200*100000*units.Pa,amount_basis=idaesHelmholtz.AmountBasis.MOLE)},
        "pressure":{0:200*100000},
    }
    seq.set_guesses_for(m.fs.main_comp.inlet,tear_guesses)

    def custom_initializer_new(unit):
        #unit.display()
        if isinstance(unit,(liqPipe,gasPipe)):
            #unit.display()
            unit.initialize(tee=False,display_after=True)
            #input('paused')
            return
        elif isinstance(unit,wellpattern_HH):
            #unit.display()
            unit.initialize(tee=False,display_after=True)
            #input('paused')
            return
        elif unit.name == 'fs.processing_facility':
            #unit.display()
            unit.initialize(display_after=False)
            #unit.display()
            #input('paused')
            return
        elif isinstance(unit,mixer):
            #unit.display()
            # ipopt.options['max_iter']=30
            # ipopt.options['acceptable_tol']=1e-2
            # set_scaling_factor(unit.enthalpy_mixing_equations,1e-3)
            # scaled_unit = pyo.TransformationFactory('core.scale_model').create_using(unit)
            # ipopt.solve(scaled_unit,tee=False)
            # pyo.TransformationFactory('core.scale_model').propagate_solution(scaled_unit,unit)
            # ipopt.options['max_iter']=3000
            # ipopt.options['acceptable_tol']=1e-6
            #unit.outlet.temperature[0] = 310
            unit.initialize()
            unit.display()
            #input('paused')
            return
        elif isinstance(unit,splitter):
            #unit.display()
            # ipopt.options['max_iter']=30
            # ipopt.options['acceptable_tol']=1e-2
            # scaled_unit = pyo.TransformationFactory('core.scale_model').create_using(unit)
            # ipopt.solve(scaled_unit,tee=False)
            # pyo.TransformationFactory('core.scale_model').propagate_solution(scaled_unit,unit)
            # ipopt.options['max_iter']=3000
            # ipopt.options['acceptable_tol']=1e-6
            unit.initialize()
            unit.display()
            #input('paused')
            return
        elif isinstance(unit,idaesPressureChanger.PressureChanger):
            unit.display()
            unit.initialize()
            unit.costing.initialize()
            #input('paused')
            ipopt.options['max_iter']=300
            ipopt.options['acceptable_tol']=1e-2
            scaled_unit = pyo.TransformationFactory('core.scale_model').create_using(unit)
            ipopt.solve(scaled_unit,tee=False)
            pyo.TransformationFactory('core.scale_model').propagate_solution(scaled_unit,unit)
            ipopt.options['max_iter']=3000
            ipopt.options['acceptable_tol']=1e-6
            #unit.display()
            #input('paused')
            return
        else:
            unit.display()
            unit.initialize()
            #input('paused')
            return

    seq.run(m,custom_initializer_new)
    
    print(f'initialization 1 done')

initialization_type1(m)
#m.fs.well1.inlet.pressure[0].fix(450*100000)
#m.fs.well2.inlet.pressure[0].fix(475*100000)

m.fs.co2_price = pyo.Param(initialize=0.068) #$/kg
m.fs.electricity_price = pyo.Param(initialize=5E-8) #$/J
#m.fs.refrigeration_price = pyo.Param(initialize=3.9E-8) #$/J
m.fs.capex_recovery_factor = pyo.Param(initialize=0.1315)

m.fs.capex = pyo.Expression(expr=
                                  m.fs.processing_facility.costing.capital_cost
                                  +m.fs.main_comp.costing.capital_cost
                                  +m.fs.recycle_comp.costing.capital_cost
                                  +m.fs.pipe0.costing.capital_cost
                                  +m.fs.pipe1.costing.capital_cost
                                  +m.fs.pipe2.costing.capital_cost
                                  +m.fs.pipe3.costing.capital_cost
                                  +m.fs.pipe4.costing.capital_cost)

m.fs.opex = pyo.Expression(expr=
                                 m.fs.electricity_price *(
                                     m.fs.main_comp.work_mechanical[0]
                                     +m.fs.recycle_comp.work_mechanical[0])*3600*24*365
                                  + m.fs.processing_facility.costing.opex)
m.fs.raw_mats = pyo.Expression(expr=m.fs.pipe0.control_volume.properties_in[0].flow_mass*m.fs.co2_price*3600*24*365)
m.fs.total_oil_production = pyo.Expression(expr=
                                           (m.fs.well1.q_OIL_PROD
                                           +m.fs.well2.q_OIL_PROD)*3600*24)

m.fs.production_target = pyo.Constraint(expr=m.fs.total_oil_production>= 450000 / 365)
m.fs.obj = pyo.Objective(expr=
                         m.fs.capex_recovery_factor*m.fs.capex
                         +m.fs.opex
                         +m.fs.raw_mats)
set_scaling_factor(m.fs.obj,1e-7)

m.fs.pipe0.inlet.flow_mol[0].unfix()
m.fs.recycle_comp.outlet.pressure[0].lb = 80 *100000
m.fs.well1.inlet.pressure[0].unfix()
m.fs.well2.inlet.pressure[0].unfix()
m.fs.main_comp.outlet.pressure[0].unfix()
m.fs.recycle_comp.outlet.pressure[0].unfix()
m.fs.purge_splitter.to_mix0.pressure[0].unfix()
m.fs.pipe0.Pdrop_constraint.deactivate()
m.fs.pipe1.Pdrop_constraint.deactivate()
m.fs.pipe2.Pdrop_constraint.deactivate()
m.fs.pipe3.Pdrop_constraint.deactivate()
m.fs.pipe4.Pdrop_constraint.deactivate()
m.fs.mix0.outlet.pressure[0].unfix()
m.fs.mix1.outlet.pressure[0].unfix()
m.fs.split1.to_well1.pressure[0].unfix()
m.fs.split1.to_pipe2.pressure[0].unfix()

with open('temps/flowsheet7_HH_postinitialization_display.txt', 'w') as f:
    with contextlib.redirect_stdout(f):
        m.display()

from co2_eor.utility.util_funcs import report_infeasibilities
report_infeasibilities(m, "temps/presolve_violations.log")

from idaes.core.scaling.autoscaling import AutoScaler
autoScaler=AutoScaler(overwrite=True)
autoScaler.scale_variables_by_magnitude(m)
autoScaler.scale_constraints_by_jacobian_norm(m)

ipopt.options['max_iter'] = 2000
ipopt.options['linear_solver']='ma27'
ipopt.options['acceptable_tol']=1E-6
ipopt.options['tol']=1E-8
ipopt.options['expect_infeasible_problem']='yes'
#ipopt.options['nlp_scaling_method']='gradient-based'
#ipopt.options['expect_infeasible_problem']='yes'
#ipopt.options['OF_evaluate_orig_obj_at_resto_trial']='no'
#ipopt.options['OF_start_with_resto']='yes'
#scale model
scaled_m = pyo.TransformationFactory("core.scale_model").create_using(m)
#solve
res=ipopt.solve(scaled_m,tee=True)
#unscale model
pyo.TransformationFactory("core.scale_model").propagate_solution(scaled_m,m)

with open('temps/flowsheet7_HH_postsolve_display.txt', 'w') as f:
    with contextlib.redirect_stdout(f):
        m.display()

report_infeasibilities(m, "temps/postsolve_violations.log")

pyo.assert_optimal_termination(res)

print(f'{pyo.value(m.fs.capex)=}')
print(f'{pyo.value(m.fs.opex)=}')
print(f'{pyo.value(m.fs.raw_mats)=}')
print(f'{pyo.value(m.fs.total_oil_production)=}')

from co2_eor.utility.util_funcs import export_flowsheet_to_excel
export_flowsheet_to_excel(m.fs, 'temps/flowsheet7_HH_solve.xlsx')