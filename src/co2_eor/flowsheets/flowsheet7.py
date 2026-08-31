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

from co2_eor.SSLW import SSLWCosting, SSLWCostingData
from co2_eor.SSLW import CompressorType, CompressorDriveType, CompressorMaterial
from co2_eor.SSLW import VesselMaterial

from co2_eor import pipeline, mixer,splitter, wellpad, heater
from co2_eor.util_funcs import export_compressor_df

from idaes.models.properties.modular_properties.base.generic_property import GenericParameterBlock
from co2_eor.liq_pipe import liqPipe
from co2_eor.wellpattern_HH import wellpattern
import co2_eor.MPF.thermo_config as thermo_config
from co2_eor.gas_pipe import gasPipe
from co2_eor.processing_facility import processingFacility
from co2_eor.mpf_to_helmholtz import mpf_helmholtz_converter

from co2_eor.gas_processing_flowsheets.absorption_plant_r2 import absorption_plant_fs

from co2_eor.util_funcs import ipopt, conopt

m = pyo.ConcreteModel()
m.fs = idaes.core.FlowsheetBlock(dynamic=False)
m.fs.props_helmholtz = idaesHelmholtz.HelmholtzParameterBlock(
        pure_component="CO2",amount_basis=idaesHelmholtz.AmountBasis.MOLE,
        state_vars=idaesHelmholtz.StateVars.PH,phase_presentation=idaesHelmholtz.PhaseType.MIX,
        #has_phase_equilibrium=False,
        )
m.fs.props_mix_vap_cubic = GenericParameterBlock(**thermo_config.configuration_vap_cubic)

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
    "GOR":4.42,
    "PC_A":0.4,
    "PC_B":0.6,
    "GB_A":0.62,
    "GB_B":1.12,
    "kovr":2.93E-14,
    "use_correction_factor":False,
    "depth":1000,
    "BHP_max":650*100000, #650?
    "outlet_pressure":30*100000,
    "outlet_temperature":300,
}

heater_config = {
    "property_package":m.fs.props_helmholtz,
}

gasPipe_config = {
    "property_package":m.fs.props_mix_vap_cubic,
    "average_pressure_type":'nonlinear',
    "allowable_stress":1380*100000,
}

#build units
m.fs.pipe0 = liqPipe(**liqPipe_config,length=15000)
#m.fs.pipe0.diameter.fix(5)
m.fs.pipe0.roughness.fix(0.0475e-3)
m.fs.pipe0.max_pressure.fix(140*100000)

m.fs.mix0 = mixer(**mixer_config,property_package=m.fs.props_helmholtz,inlet_list=['from_pipe0','from_purge_splitter'])

m.fs.main_comp = idaesPressureChanger.PressureChanger(**compressor_config)
m.fs.main_comp.efficiency_isentropic.fix(0.85)

m.fs.chiller = heater(**heater_config)

m.fs.pipe1 = liqPipe(**liqPipe_config,length=3000)
#m.fs.pipe1.diameter.fix(5)
m.fs.pipe1.roughness.fix(0.0475e-3)
m.fs.pipe1.max_pressure.fix(650*100000)

m.fs.split1 = splitter(**splitter_config,outlet_list=['to_well1','to_pipe2'])

m.fs.well1 = wellpattern(**wellpattern_config)
m.fs.well1.HCPV.fix(1.1)

m.fs.pipe2 = liqPipe(**liqPipe_config,length=2200)
#m.fs.pipe2.diameter.fix(5)
m.fs.pipe2.roughness.fix(0.0475e-3)
m.fs.pipe2.max_pressure.fix(650*100000)

m.fs.well2 = wellpattern(**wellpattern_config)
m.fs.well2.HCPV.fix(0.9)

m.fs.pipe3 = gasPipe(**gasPipe_config,length=2300)
#m.fs.pipe3.diameter.fix(5)
m.fs.pipe3.roughness.fix(0.0475e-3)
m.fs.pipe3.max_pressure.fix(50*100000)

m.fs.mix1 = mixer(**mixer_config,property_package=m.fs.props_mix_vap_cubic,inlet_list=['from_well1','from_pipe3'])

m.fs.pipe4 = gasPipe(**gasPipe_config,length=3100)
#m.fs.pipe4.diameter.fix(5)
m.fs.pipe4.roughness.fix(0.0473e-3)
m.fs.pipe4.max_pressure.fix(50*100000)
"""
m.fs.processing_facility = absorption_plant_fs(
    property_package=m.fs.props_mix_vap_cubic,
    solutes=['co2','ch4'],
    henry_coefficients={
        ('co2','A'):13.828+6.9,('co2','B'):-1720,
        ('ch4','A'):16.531+6.9,('ch4','B'):-1720,
    },
)
m.fs.processing_facility.column.num_trays.fix(10)
m.fs.processing_facility.column.top_inlet.flow_mol_comp[0,'co2'].fix(0)
m.fs.processing_facility.column.top_inlet.flow_mol_comp[0,'ch4'].fix(0)
m.fs.processing_facility.column.top_inlet.flow_mol_comp[0,'selexol'].fix(1000)
m.fs.processing_facility.column.top_inlet.pressure[0].fix(2*100000)
m.fs.processing_facility.column.top_inlet.temperature[0].fix(310)
"""
m.fs.processing_facility=processingFacility(
    property_package=m.fs.props_mix_vap_cubic,
    key_component='co2',
    minimum_inlet_pressure=5*100000,
    recycle_pressure=5*100000,
)

m.fs.eos_converter = mpf_helmholtz_converter(
    property_package_in = m.fs.props_mix_vap_cubic,
    property_package_out = m.fs.props_helmholtz,
    conversion_type = 'convert_all_mol'
)

m.fs.recycle_comp = idaesPressureChanger.PressureChanger(**compressor_config)
m.fs.recycle_comp.efficiency_isentropic.fix(0.85)

m.fs.purge_splitter = splitter(**splitter_config,outlet_list=['to_purge','to_mix0'])
m.fs.purge_splitter.positive_purge_eq = pyo.Constraint(expr=m.fs.purge_splitter.to_purge.flow_mol[0]>=0) #enforce material does not go into system from purge stream
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
#m.fs.s_pipe4_processing_facility = Arc(source=m.fs.pipe4.outlet,destination=m.fs.processing_facility.translator1.inlet)
#m.fs.s_processing_facility_eos_converter = Arc(source=m.fs.processing_facility.translator2.outlet,destination=m.fs.eos_converter.inlet)
m.fs.s_pipe4_processing_facility = Arc(source=m.fs.pipe4.outlet,destination=m.fs.processing_facility.inlet)
m.fs.s_processing_facility_eos_converter = Arc(source=m.fs.processing_facility.recycle,destination=m.fs.eos_converter.inlet)
m.fs.s_eos_converter_recycle_comp = Arc(source=m.fs.eos_converter.outlet,destination=m.fs.recycle_comp.inlet)
m.fs.s_recycle_comp_purge_splitter = Arc(source=m.fs.recycle_comp.outlet,destination=m.fs.purge_splitter.inlet)
m.fs.s_purge_splitter_mix0 = Arc(source=m.fs.purge_splitter.to_mix0,destination=m.fs.mix0.from_purge_splitter)

pyo.TransformationFactory("network.expand_arcs").apply_to(m)

#fix degrees of freedom
m.fs.pipe0.inlet.pressure[0].fix(100*100000)
m.fs.pipe0.inlet.enth_mol[0].fix(m.fs.props_helmholtz.htpx(T=310*units.K,p=100*100000*units.Pa,amount_basis=idaesHelmholtz.AmountBasis.MOLE))

m.fs.main_comp.outlet.pressure[0].fix(500*100000)
m.fs.chiller.outlet_temp_eq = pyo.Constraint(expr=m.fs.chiller.control_volume.properties_out[0].temperature==310)

m.fs.recycle_comp.outlet.pressure[0].fix(86*100000)

m.fs.purge_splitter.split_fraction[0,'to_purge'].fix(1e-14)
m.fs.purge_splitter.to_mix0.pressure[0].fix(85*100000)

m.fs.well1.inlet.pressure[0].fix(450*100000)
m.fs.well2.inlet.pressure[0].fix(475*100000)

#fix pressure drops instead of diameters
m.fs.pipe0.Pdrop_constraint = pyo.Constraint(expr=m.fs.pipe0.Pdrop==15*100000)
m.fs.pipe1.Pdrop_constraint = pyo.Constraint(expr=m.fs.pipe1.Pdrop==15*100000)
m.fs.pipe2.Pdrop_constraint = pyo.Constraint(expr=m.fs.pipe2.Pdrop==8*100000)
m.fs.pipe3.Pdrop_constraint = pyo.Constraint(expr=m.fs.pipe3.Pdrop==16*100000)
m.fs.pipe4.Pdrop_constraint = pyo.Constraint(expr=m.fs.pipe4.Pdrop==6*100000)

m.fs.mix0.outlet.pressure[0].fix(84*100000)
m.fs.mix1.outlet.pressure[0].fix(13.5*100000)

#ancillary degrees of freedom
#m.fs.split1.to_well1.pressure[0].fix(450*100000)
#m.fs.split1.to_pipe2.pressure[0].fix(530*100000)
#m.fs.mix0.outlet.pressure[0].fix(79*100000)
#m.fs.mix1.outlet.pressure[0].fix(8*100000)

#check degrees of freedom
print(f'number of variables={len(list(m.component_data_objects(pyo.Var)))}')
print(f'number of constraints={len(list(m.component_data_objects(pyo.Constraint)))}')
DoF = idaes.core.util.model_statistics.degrees_of_freedom(m)
print(f'degrees of freedom={DoF}')
with open('temps/flowsheet_7_preinitialization_pprint.txt', 'w') as f:
    with contextlib.redirect_stdout(f):
        m.pprint()
print("Press any key to continue (or 'q' to quit)...")
key = readchar.readkey()
if key.lower() == "q":
    print("Exiting program.")
    sys.exit()
print(f"Resumed after pressing: {key}")

"""
m.fs.well1.inlet.pressure[0].fix(450*100000)
m.fs.well1.inlet.enth_mol[0].fix(m.fs.props_helmholtz.htpx(T=310*units.K,p=450*100000*units.Pa,amount_basis=idaesHelmholtz.AmountBasis.MOLE))
m.fs.well1.initialize(display_after=True)
input('paused')

m.fs.well2.inlet.pressure[0].fix(550*100000)
m.fs.well2.inlet.enth_mol[0].fix(m.fs.props_helmholtz.htpx(T=310*units.K,p=550*100000*units.Pa,amount_basis=idaesHelmholtz.AmountBasis.MOLE))
m.fs.well2.initialize(display_after=True)
input('paused')
"""

#initialization degrees of freedom
m.fs.pipe0.inlet.flow_mol[0].fix(40)
m.fs.well1.inlet.pressure[0].unfix()
m.fs.well2.inlet.pressure[0].unfix()

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

    tear_guesses = {
        "flow_mol":{
            (0):60,
        },
        "enth_mol":{0:m.fs.props_helmholtz.htpx(T=300*units.K,p=80*100000*units.Pa,amount_basis=idaesHelmholtz.AmountBasis.MOLE)},
        "pressure":{0:80*100000},
    }
    seq.set_guesses_for(m.fs.mix0.outlet,tear_guesses)

    from co2_eor.util_funcs import custom_initialize_unit
    def custom_initializer_new(unit):
        #unit.display()
        if isinstance(unit,(liqPipe,gasPipe)):
            #unit.display()
            unit.initialize(tee=False,display_after=False)
            return
        elif isinstance(unit,wellpattern):
            unit.initialize(display_after=False)
            return
        elif unit.name == 'fs.processing_facility':
            #unit.display()
            unit.initialize(display_after=False)
            return
        elif isinstance(unit,mixer):
            #unit.display()
            ipopt.options['max_iter']=30
            ipopt.options['acceptable_tol']=1e-2
            set_scaling_factor(unit.enthalpy_mixing_equations,1e-3)
            scaled_unit = pyo.TransformationFactory('core.scale_model').create_using(unit)
            ipopt.solve(scaled_unit,tee=False)
            pyo.TransformationFactory('core.scale_model').propagate_solution(scaled_unit,unit)
            ipopt.options['max_iter']=3000
            ipopt.options['acceptable_tol']=1e-6
            return
        elif isinstance(unit,idaesPressureChanger.PressureChanger):
            unit.initialize()
            return
        else:
            #unit.display()
            unit.initialize()
            return

    seq.run(m,custom_initializer_new)
    
    print(f'initialization 1 done')

initialization_type1(m)

m.fs.pipe0.inlet.flow_mol[0].unfix()
m.fs.well1.inlet.pressure[0].fix(450*100000)
m.fs.well2.inlet.pressure[0].fix(475*100000)


def initialization_2(m):
    m.fs.well1.inlet.pressure[0].fix(450*100000)
    m.fs.well1.inlet.enth_mol[0].fix(m.fs.props_helmholtz.htpx(T=310*units.K,p=450*100000*units.Pa,amount_basis=idaesHelmholtz.AmountBasis.MOLE))

    m.fs.well2.inlet.pressure[0].fix(550*100000)
    m.fs.well2.inlet.enth_mol[0].fix(m.fs.props_helmholtz.htpx(T=310*units.K,p=550*100000*units.Pa,amount_basis=idaesHelmholtz.AmountBasis.MOLE))
    
    from pyomo.network import SequentialDecomposition
    seq = SequentialDecomposition()
    G = seq.create_graph(m)
    custom_order = [
        [m.fs.well1],
        [m.fs.well2],
        [m.fs.pipe3],
        [m.fs.mix1],
        [m.fs.pipe4],
        [m.fs.processing_facility],
        [m.fs.eos_converter],
        [m.fs.recycle_comp],
        [m.fs.purge_splitter]
    ]
    def custom_initializer_new(unit):
        if isinstance(unit,(liqPipe,gasPipe)):
            #unit.display()
            unit.initialize(tee=False)
            return
        if isinstance(unit,wellpattern):
            unit.initialize(display_after=False)
            return
        if unit.name == 'fs.processing_facility.column':
            unit.initialize()
            return
        if isinstance(unit,mixer):
            ipopt.options['max_iter']=10
            ipopt.solve(unit)
            ipopt.options['max_iter']=3000
        else:
            unit.initialize()
            return

    seq.run_order(G,order=custom_order,function=custom_initializer_new)

#initialization_2(m)


with open('temps/flowsheet7_postinitialization_display.txt', 'w') as f:
    with contextlib.redirect_stdout(f):
        m.display()

ipopt.options['linear_solver']='ma27'
#scale model
scaled_m = pyo.TransformationFactory("core.scale_model").create_using(m)
#solve
res=ipopt.solve(scaled_m,tee=True)
#unscale model
pyo.TransformationFactory("core.scale_model").propagate_solution(scaled_m,m)

with open('temps/flowsheet7_postsolve_display.txt', 'w') as f:
    with contextlib.redirect_stdout(f):
        m.display()

pyo.assert_optimal_termination(res)