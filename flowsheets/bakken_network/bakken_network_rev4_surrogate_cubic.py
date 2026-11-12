import pyomo.environ as pyo
from pyomo.environ import units as pyunits
import idaes.core
from pyomo.network import Arc
import contextlib
import pickle
import time
from pathlib import Path

import idaes.models.properties.general_helmholtz as idaes_helmholtz
import idaes.models.unit_models.pressure_changer as idaes_pressure_changer

from idaes.core.util.scaling import set_scaling_factor

from co2_eor import mixer, splitter, heater
from idaes.models.properties.modular_properties.base.generic_property import GenericParameterBlock
from co2_eor.liq_pipe import liqPipe
from co2_eor.wellpattern_unit import wellpattern
import co2_eor.MPF.thermo_config as thermo_config
from co2_eor.gas_pipe import gasPipe

from co2_eor.surrogates.rectisol_surrogate import rectisolSurrogate
from co2_eor.mpf_to_helmholtz import mpf_helmholtz_converter

from co2_eor.SSLW import SSLWCosting, SSLWCostingData

from co2_eor.utility.solvers import ipopt

m_name = 'bakken_network_rev4_surrogate_cubic'
output_dir = 'temps/' + m_name + '/'
Path(output_dir).mkdir(parents=True, exist_ok=True)
m = pyo.ConcreteModel()
m.fs = idaes.core.FlowsheetBlock(dynamic=False)

m.fs.props_helmholtz = idaes_helmholtz.HelmholtzParameterBlock(
        pure_component="CO2",amount_basis=idaes_helmholtz.AmountBasis.MOLE,
        state_vars=idaes_helmholtz.StateVars.PH,phase_presentation=idaes_helmholtz.PhaseType.MIX,
        #has_phase_equilibrium=False,
        )
m.fs.props_mix_vap_cubic = GenericParameterBlock(**thermo_config.configuration_vap_cubic)
m.fs.props_mix_liq_cubic = GenericParameterBlock(**thermo_config.configuration_liq_cubic)

m.fs.costing = SSLWCosting()

#define default unit configurations
compressor_config = {
    "property_package":m.fs.props_mix_liq_cubic,
    "dynamic":False,
    "compressor":True,
    "thermodynamic_assumption":idaes_pressure_changer.ThermodynamicAssumption.isentropic
}

liqpipe_config = {
        "property_package":m.fs.props_mix_liq_cubic,
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
        "momentum_mixing_type":4 #for inequality mixing type
}

splitter_config = {
        "property_package":m.fs.props_mix_liq_cubic,
        "ideal_separation":False,
        "energy_split_basis":2, #for equal_molar_enthalpy
        "momentum_balance_type":0, #for none
}

wellpattern_config = {
    "primary_property_package":m.fs.props_mix_liq_cubic,
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
    "multiplier":1, 
}

heater_config = {
    "property_package":m.fs.props_mix_liq_cubic,
}

gaspipe_config = {
    "property_package":m.fs.props_mix_vap_cubic,
    "average_pressure_type":'nonlinear',
    "allowable_stress":1380*100000,
}

liqpipe_costing_config = {
    "flowsheet_costing_block":m.fs.costing,
    "costing_method":SSLWCostingData.cost_horizontal_vessel,
    "costing_method_arguments":{
        "material_type":'StainlessSteel316',
        "include_platforms_ladders":False,
    }
}

gaspipe_costing_config = {
    "flowsheet_costing_block":m.fs.costing,
    "costing_method":SSLWCostingData.cost_horizontal_vessel,
    "costing_method_arguments":{
        "material_type":'StainlessSteel316',
        "include_platforms_ladders":False,
    }
}

compressor_costing_config={
    "flowsheet_costing_block":m.fs.costing,
    "costing_method":SSLWCostingData.cost_compressor,
    "costing_method_arguments":{
        "compressor_type":'Centrifugal',
        "drive_type":'ElectricMotor',
        "material_type":'StainlessSteel',
    }
}

surrogate_dir = 'temps/rectisol_surrogate/'
surrogate_object = idaes.core.surrogate.AlamoSurrogate.load_from_file(surrogate_dir+'alamo_surrogate.json')
with open(surrogate_dir+"input_bounds.pkl", "rb") as file:
    surrogate_input_bounds = pickle.load(file)
with open(surrogate_dir+"hyperplane_list.pkl", "rb") as file:
    hyperplane_list = pickle.load(file)
with open(surrogate_dir+"co2_recovery_min.pkl", "rb") as file:
    co2_recovery_min = pickle.load(file)
processingfacility_config = {
    "inlet_property_package":m.fs.props_mix_vap_cubic,
    "outlet_property_package":m.fs.props_mix_vap_cubic,
    "recycle_property_package":m.fs.props_mix_liq_cubic,
    "surrogate_input_bounds":surrogate_input_bounds,
    "unit_model_parameters":{
        "column":{
            "stage_height":0.6067,
            "allowable_stress":1300e5,
        },
        "flash3":{
            "residence_time":60,
            "allowable_stress":1300e5,
        },
        "flash2":{
            "residence_time":60,
            "allowable_stress":1300e5,
        },
        "flash1":{
            "residence_time":60,
            "allowable_stress":1300e5,
        },
        "solvent_tank":{
            "storage_time":60*3600*24,
        },
    }, 
    "global_parameters":{
        "solvent_MW":32.04e-3,
        "solvent_mass_dens":800,
        "recycle_outlet_temperature":350,
        "co2_recovery_min":co2_recovery_min
    },
    "surrogate_object":surrogate_object,
    "hyperplane_cuts":hyperplane_list,
}

#distribution network
#source branch
m.fs.source_comp = idaes_pressure_changer.PressureChanger(**compressor_config)
m.fs.source_comp.costing = idaes.core.UnitModelCostingBlock(**compressor_costing_config)
m.fs.source_comp.efficiency_isentropic.fix(0.85)

m.fs.pipe0 = liqPipe(**liqpipe_config,length=15000)
m.fs.pipe0.costing = idaes.core.UnitModelCostingBlock(**liqpipe_costing_config)
m.fs.pipe0.diameter = 1.2
m.fs.pipe0.roughness.fix(0.0475e-3)

m.fs.recycle_comp = idaes_pressure_changer.PressureChanger(**compressor_config)
m.fs.recycle_comp.costing = idaes.core.UnitModelCostingBlock(**compressor_costing_config)
m.fs.recycle_comp.efficiency_isentropic.fix(0.85)

m.fs.chiller = heater(**heater_config)

m.fs.s_source_comp_pipe0 = Arc(source=m.fs.source_comp.outlet,destination=m.fs.pipe0.inlet)

m.fs.mix0 = mixer(**mixer_config,property_package=m.fs.props_mix_liq_cubic,inlet_list=['from_pipe0','from_recycle_comp'])

m.fs.s_recycle_comp_mix0 = Arc(source=m.fs.recycle_comp.outlet,destination=m.fs.mix0.from_recycle_comp)
m.fs.s_pipe0_mix0 = Arc(source=m.fs.pipe0.outlet,destination=m.fs.mix0.from_pipe0)
m.fs.s_mix0_chiller = Arc(source=m.fs.mix0.outlet,destination=m.fs.chiller.inlet)

m.fs.split0 = splitter(**splitter_config,outlet_list=['to_pipe1','to_pipe3','to_pipe6'])

m.fs.s_chiller_split0 = Arc(source=m.fs.chiller.outlet,destination=m.fs.split0.inlet)

#first branch, wells 1 and 2
m.fs.pipe1 = liqPipe(**liqpipe_config,length=800)
m.fs.pipe1.costing = idaes.core.UnitModelCostingBlock(**liqpipe_costing_config)
m.fs.pipe1.diameter = 1
m.fs.pipe1.roughness.fix(0.0475e-3)

m.fs.split1 = splitter(**splitter_config,outlet_list=['to_comp1','to_pipe2'])

m.fs.comp1 = idaes_pressure_changer.PressureChanger(**compressor_config)
m.fs.comp1.costing = idaes.core.UnitModelCostingBlock(**compressor_costing_config)
m.fs.comp1.efficiency_isentropic.fix(0.85)
m.fs.comp1.pressure_minimum_eq = pyo.Constraint(expr=m.fs.comp1.inlet.pressure[0]>=80*100000)

m.fs.well1 = wellpattern(**wellpattern_config)
m.fs.well1.HCPV.fix(1.1)

m.fs.pipe2 = liqPipe(**liqpipe_config,length=1600)
m.fs.pipe2.costing = idaes.core.UnitModelCostingBlock(**liqpipe_costing_config)
m.fs.pipe2.diameter = 1
m.fs.pipe2.roughness.fix(0.0475e-3)

m.fs.comp2 = idaes_pressure_changer.PressureChanger(**compressor_config)
m.fs.comp2.costing = idaes.core.UnitModelCostingBlock(**compressor_costing_config)
m.fs.comp2.efficiency_isentropic.fix(0.85)
m.fs.comp2.pressure_minimum_eq = pyo.Constraint(expr=m.fs.comp2.inlet.pressure[0]>=80*100000)

m.fs.well2 = wellpattern(**wellpattern_config)
m.fs.well2.HCPV.fix(0.9)

m.fs.s_split0_pipe1 = Arc(source=m.fs.split0.to_pipe1,destination=m.fs.pipe1.inlet)
m.fs.s_pipe1_split1 = Arc(source=m.fs.pipe1.outlet,destination=m.fs.split1.inlet)
m.fs.s_split1_comp1 = Arc(source=m.fs.split1.to_comp1,destination=m.fs.comp1.inlet)
m.fs.s_comp1_well1 = Arc(source=m.fs.comp1.outlet,destination=m.fs.well1.inlet)
m.fs.s_split1_pipe2 = Arc(source=m.fs.split1.to_pipe2,destination=m.fs.pipe2.inlet)
m.fs.s_pipe2_comp2 = Arc(source=m.fs.pipe2.outlet,destination=m.fs.comp2.inlet)
m.fs.s_comp2_well2 = Arc(source=m.fs.comp2.outlet,destination=m.fs.well2.inlet)

#second branch, wells 3 and 4
m.fs.pipe3 = liqPipe(**liqpipe_config,length=3200)
m.fs.pipe3.costing=idaes.core.UnitModelCostingBlock(**liqpipe_costing_config)
m.fs.pipe3.diameter = 1
m.fs.pipe3.roughness.fix(0.0475e-3)

m.fs.split2 = splitter(**splitter_config,outlet_list=['to_pipe4','to_pipe5'])

m.fs.pipe4 = liqPipe(**liqpipe_config,length=4000)
m.fs.pipe4.costing=idaes.core.UnitModelCostingBlock(**liqpipe_costing_config)
m.fs.pipe4.diameter = 1
m.fs.pipe4.roughness.fix(0.0475e-3)

m.fs.comp3 = idaes_pressure_changer.PressureChanger(**compressor_config)
m.fs.comp3.costing = idaes.core.UnitModelCostingBlock(**compressor_costing_config)
m.fs.comp3.efficiency_isentropic.fix(0.85)
m.fs.comp3.pressure_minimum_eq = pyo.Constraint(expr=m.fs.comp3.inlet.pressure[0]>=80*100000)

m.fs.well3 = wellpattern(**wellpattern_config)
m.fs.well3.HCPV.fix(0.8)

m.fs.pipe5 = liqPipe(**liqpipe_config,length=2400)
m.fs.pipe5.costing=idaes.core.UnitModelCostingBlock(**liqpipe_costing_config)
m.fs.pipe5.diameter = 1
m.fs.pipe5.roughness.fix(0.0475e-3)

m.fs.comp4 = idaes_pressure_changer.PressureChanger(**compressor_config)
m.fs.comp4.costing = idaes.core.UnitModelCostingBlock(**compressor_costing_config)
m.fs.comp4.efficiency_isentropic.fix(0.85)
m.fs.comp4.pressure_minimum_eq = pyo.Constraint(expr=m.fs.comp4.inlet.pressure[0]>=80*100000)

m.fs.well4 = wellpattern(**wellpattern_config)
m.fs.well4.HCPV.fix(0.7)

m.fs.s_split0_pipe3 = Arc(source=m.fs.split0.to_pipe3,destination=m.fs.pipe3.inlet)
m.fs.s_pipe3_split2 = Arc(source=m.fs.pipe3.outlet,destination=m.fs.split2.inlet)
m.fs.s_split2_pipe4 = Arc(source=m.fs.split2.to_pipe4,destination=m.fs.pipe4.inlet)
m.fs.s_pipe4_comp3 = Arc(source=m.fs.pipe4.outlet,destination=m.fs.comp3.inlet)
m.fs.s_comp3_well3 = Arc(source=m.fs.comp3.outlet,destination=m.fs.well3.inlet)
m.fs.s_split2_pipe5 = Arc(source=m.fs.split2.to_pipe5,destination=m.fs.pipe5.inlet)
m.fs.s_pipe5_comp4 = Arc(source=m.fs.pipe5.outlet,destination=m.fs.comp4.inlet)
m.fs.s_comp4_well4 = Arc(source=m.fs.comp4.outlet,destination=m.fs.well4.inlet)

#third branch, wells 5 and 6
m.fs.pipe6 = liqPipe(**liqpipe_config,length=3200)
m.fs.pipe6.costing=idaes.core.UnitModelCostingBlock(**liqpipe_costing_config)
m.fs.pipe6.diameter = 1
m.fs.pipe6.roughness.fix(0.0475e-3)

m.fs.split3 = splitter(**splitter_config,outlet_list=['to_pipe7','to_pipe8'])

m.fs.pipe7 = liqPipe(**liqpipe_config,length=6900)
m.fs.pipe7.costing=idaes.core.UnitModelCostingBlock(**liqpipe_costing_config)
m.fs.pipe7.diameter = 1
m.fs.pipe7.roughness.fix(0.0475e-3)

m.fs.comp5 = idaes_pressure_changer.PressureChanger(**compressor_config)
m.fs.comp5.costing = idaes.core.UnitModelCostingBlock(**compressor_costing_config)
m.fs.comp5.efficiency_isentropic.fix(0.85)
m.fs.comp5.pressure_minimum_eq = pyo.Constraint(expr=m.fs.comp5.inlet.pressure[0]>=80*100000)

m.fs.well5 = wellpattern(**wellpattern_config)
m.fs.well5.HCPV.fix(1.0)

m.fs.pipe8 = liqPipe(**liqpipe_config,length=4000)
m.fs.pipe8.costing=idaes.core.UnitModelCostingBlock(**liqpipe_costing_config)
m.fs.pipe8.diameter = 1
m.fs.pipe8.roughness.fix(0.0475e-3)

m.fs.comp6 = idaes_pressure_changer.PressureChanger(**compressor_config)
m.fs.comp6.costing = idaes.core.UnitModelCostingBlock(**compressor_costing_config)
m.fs.comp6.efficiency_isentropic.fix(0.85)
m.fs.comp6.pressure_minimum_eq = pyo.Constraint(expr=m.fs.comp6.inlet.pressure[0]>=80*100000)

m.fs.well6 = wellpattern(**wellpattern_config)
m.fs.well6.HCPV.fix(0.6)

m.fs.s_split0_pipe6 = Arc(source=m.fs.split0.to_pipe6,destination=m.fs.pipe6.inlet)
m.fs.s_pipe6_split3 = Arc(source=m.fs.pipe6.outlet,destination=m.fs.split3.inlet)
m.fs.s_split3_pipe7 = Arc(source=m.fs.split3.to_pipe7,destination=m.fs.pipe7.inlet)
m.fs.s_pipe7_comp5 = Arc(source=m.fs.pipe7.outlet,destination=m.fs.comp5.inlet)
m.fs.s_comp5_well5 = Arc(source=m.fs.comp5.outlet,destination=m.fs.well5.inlet)
m.fs.s_split3_pipe8 = Arc(source=m.fs.split3.to_pipe8,destination=m.fs.pipe8.inlet)
m.fs.s_pipe8_comp6 = Arc(source=m.fs.pipe8.outlet,destination=m.fs.comp6.inlet)
m.fs.s_comp6_well6 = Arc(source=m.fs.comp6.outlet,destination=m.fs.well6.inlet)

#collection network
m.fs.mixPF = mixer(**mixer_config,property_package=m.fs.props_mix_vap_cubic,inlet_list=['from_pipe9','from_pipe12','from_pipe14'])

#collection from wells 1,2, 4
m.fs.pipe9 = gasPipe(**gaspipe_config,length=1000)
m.fs.pipe9.costing=idaes.core.UnitModelCostingBlock(**gaspipe_costing_config)
m.fs.pipe9.diameter = 1
m.fs.pipe9.roughness.fix(0.0475e-3)

m.fs.mix1 = mixer(**mixer_config,property_package=m.fs.props_mix_vap_cubic,inlet_list=['from_well1','from_pipe10'])

m.fs.pipe10 = gasPipe(**gaspipe_config,length=1600)
m.fs.pipe10.costing=idaes.core.UnitModelCostingBlock(**gaspipe_costing_config)
m.fs.pipe10.diameter = 1
m.fs.pipe10.roughness.fix(0.0475e-3)

m.fs.mix2 = mixer(**mixer_config,property_package=m.fs.props_mix_vap_cubic,inlet_list=['from_well2','from_pipe11'])

m.fs.pipe11 = gasPipe(**gaspipe_config,length=3200)
m.fs.pipe11.costing=idaes.core.UnitModelCostingBlock(**gaspipe_costing_config)
m.fs.pipe11.diameter = 1
m.fs.pipe11.roughness.fix(0.0475e-3)

m.fs.s_well4_pipe11 = Arc(source=m.fs.well4.outlet,destination=m.fs.pipe11.inlet)
m.fs.s_pipe11_mix2 = Arc(source=m.fs.pipe11.outlet,destination=m.fs.mix2.from_pipe11)
m.fs.s_well2_mix2 = Arc(source=m.fs.well2.outlet,destination=m.fs.mix2.from_well2)
m.fs.s_mix2_pipe10 = Arc(source=m.fs.mix2.outlet,destination=m.fs.pipe10.inlet)
m.fs.s_pipe10_mix1 = Arc(source=m.fs.pipe10.outlet,destination=m.fs.mix1.from_pipe10)
m.fs.s_well1_mix1 = Arc(source=m.fs.well1.outlet,destination=m.fs.mix1.from_well1)
m.fs.s_mix1_pipe9 = Arc(source=m.fs.mix1.outlet,destination=m.fs.pipe9.inlet)
m.fs.s_pipe9_mixPF = Arc(source=m.fs.pipe9.outlet,destination=m.fs.mixPF.from_pipe9)

#collection from wells 5,6
m.fs.pipe12 = gasPipe(**gaspipe_config,length=2000)
m.fs.pipe12.costing=idaes.core.UnitModelCostingBlock(**gaspipe_costing_config)
m.fs.pipe12.diameter = 1
m.fs.pipe12.roughness.fix(0.0475e-3)

m.fs.mix3 = mixer(**mixer_config,property_package=m.fs.props_mix_vap_cubic,inlet_list=['from_well6','from_pipe13'])

m.fs.pipe13 = gasPipe(**gaspipe_config,length=5500)
m.fs.pipe13.costing=idaes.core.UnitModelCostingBlock(**gaspipe_costing_config)
m.fs.pipe13.diameter = 1
m.fs.pipe13.roughness.fix(0.0475e-3)

m.fs.s_well5_pipe13 = Arc(source=m.fs.well5.outlet,destination=m.fs.pipe13.inlet)
m.fs.s_pipe13_mix3 = Arc(source=m.fs.pipe13.outlet,destination=m.fs.mix3.from_pipe13)
m.fs.s_well6_mix3 = Arc(source=m.fs.well6.outlet,destination=m.fs.mix3.from_well6)
m.fs.s_mix3_pipe12 = Arc(source=m.fs.mix3.outlet,destination=m.fs.pipe12.inlet)
m.fs.s_pipe12_mixPF = Arc(source=m.fs.pipe12.outlet,destination=m.fs.mixPF.from_pipe12)

#collection from well 3
m.fs.pipe14 = gasPipe(**gaspipe_config,length=2800)
m.fs.pipe14.costing=idaes.core.UnitModelCostingBlock(**gaspipe_costing_config)
m.fs.pipe14.diameter = 1
m.fs.pipe14.roughness.fix(0.0475e-3)

m.fs.s_well3_pipe14 = Arc(source=m.fs.well3.outlet,destination=m.fs.pipe14.inlet)
m.fs.s_pipe14_mixPF = Arc(source=m.fs.pipe14.outlet,destination=m.fs.mixPF.from_pipe14)

#processing facility to recycle compressor
m.fs.processing_facility = rectisolSurrogate(**processingfacility_config)

m.fs.purge_splitter = splitter(**splitter_config,outlet_list=['to_purge','to_recycle_comp'])
m.fs.purge_splitter.to_purge.flow_mol_comp[0,'co2'].lb = 0
m.fs.purge_splitter.to_purge.flow_mol_comp[0,'ch4'].lb = 0

m.fs.s_mixPF_processing_facility = Arc(source=m.fs.mixPF.outlet,destination=m.fs.processing_facility.inlet)
m.fs.s_processing_facility_purge_splitter = Arc(source=m.fs.processing_facility.recycle,destination=m.fs.purge_splitter.inlet)
m.fs.s_purge_splitter_recycle_comp = Arc(source=m.fs.purge_splitter.to_recycle_comp,destination=m.fs.recycle_comp.inlet)

pyo.TransformationFactory("network.expand_arcs").apply_to(m)

#network specifications
m.fs.source_comp.inlet.pressure[0].fix(100*100000)
m.fs.source_comp.inlet.temperature[0].fix(300)
m.fs.source_comp.inlet.flow_mol_comp[0,'ch4'].fix(0)

m.fs.chiller.outlet_temp_eq = pyo.Constraint(expr=m.fs.chiller.control_volume.properties_out[0].temperature==318)

with open(output_dir+m_name+'_postbuild_pprint.txt','w') as f:
    with contextlib.redirect_stdout(f):
        m.pprint()

#check degrees of freedom
DoF = idaes.core.util.model_statistics.degrees_of_freedom(m)
print(f'degrees of freedom={DoF}')

def initialization1():
    #set initialization dof
    m.fs.source_comp.outlet.pressure[0].fix(201*100000)
    m.fs.recycle_comp.outlet.pressure[0].fix(200*100000)
    m.fs.purge_splitter.split_fraction[0,'to_purge'].fix(1e-14)
    m.fs.purge_splitter.to_purge.pressure[0].fix(70*100000)
    #m.fs.purge_splitter.to_recycle_comp.pressure[0].fix(70*100000)
    #fix pressure drops instead of diameters
    for i in range(0,15):
        pipe = getattr(m.fs,"pipe"+str(i))
        if isinstance(pipe,liqPipe):
            temp = 1
            pipe.max_pressure.fix(651*100000)
        else:
            temp = 4
            pipe.max_pressure.fix(31*100000)
        pipe.Pdrop_constraint = pyo.Constraint(expr=pipe.Pdrop==temp*100000)
        set_scaling_factor(pipe.Pdrop_constraint,1e-5)
        
    #splitter outlet pressures
    for j in m.fs.split0.create_outlet_list():
        getattr(m.fs.split0,j).pressure[0].fix(195*100000)
    for j in m.fs.split1.create_outlet_list():
        getattr(m.fs.split1,j).pressure[0].fix(190*100000)
    for j in m.fs.split2.create_outlet_list():
        getattr(m.fs.split2,j).pressure[0].fix(185*100000)
    for j in m.fs.split3.create_outlet_list():
        getattr(m.fs.split3,j).pressure[0].fix(185*100000)

    #mixer outlet pressures
    m.fs.mix0.outlet.pressure[0].fix(200*100000)
    m.fs.mix1.outlet.pressure[0].fix(18*100000)
    m.fs.mix2.outlet.pressure[0].fix(22*100000)
    m.fs.mix3.outlet.pressure[0].fix(22*100000)
    m.fs.mixPF.outlet.pressure[0].fix(14*100000)

    #other compressors
    for i in range(1,7):
        getattr(m.fs, "comp"+str(i)).outlet.pressure[0].fix(650*100000)

    #processing facility d.o.f
    m.fs.processing_facility.loop_solvent_flow.fix(800)
    m.fs.processing_facility.exchanger1.outlet_temperature.fix(240)
    m.fs.processing_facility.flash3.pressure.fix(2*100000)
    m.fs.processing_facility.column.Nstages.fix(4)
    m.fs.processing_facility.deactivate_surrogate_input_bounds()

    #guess inlet flowrate
    m.fs.source_comp.inlet.flow_mol_comp[0,'co2'].fix(220)

    print(f'degrees of freedom for initialization ={idaes.core.util.model_statistics.degrees_of_freedom(m)}')

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
        "flow_mol":{(0,'co2'):500,(0,'ch4'):0},
        "enth_mol":{0:m.fs.props_helmholtz.htpx(T=300*pyunits.K,p=200*100000*pyunits.Pa,amount_basis=idaes_helmholtz.AmountBasis.MOLE)},
        "pressure":{0:200*100000},
    }
    seq.set_guesses_for(m.fs.mix0.from_recycle_comp,tear_guesses)

    def custom_initializer_func(unit):
        unit.display()
        if isinstance(unit,(liqPipe,gasPipe)):
            unit.initialize(tee=False,display_after=True)
            return
        elif isinstance(unit,wellpattern):
            unit.initialize(tee=False,display_after=True)
            # number = unit.name[-1]
            # comp = getattr(m.fs,"comp"+number)
            # unit.inlet.pressure[0].unfix()
            # unit.initialize()
            # comp.outlet.pressure[0].fix(pyo.value(unit.inlet.pressure[0]))
            # comp.initialize()
            # comp.costing.initialize()
            return
        elif isinstance(unit,rectisolSurrogate):
            unit.initialize(display_after=True)
            return
        elif isinstance(unit,mixer):
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
            # ipopt.options['max_iter']=30
            # ipopt.options['acceptable_tol']=1e-2
            # scaled_unit = pyo.TransformationFactory('core.scale_model').create_using(unit)
            # ipopt.solve(scaled_unit,tee=False)
            # pyo.TransformationFactory('core.scale_model').propagate_solution(scaled_unit,unit)
            # ipopt.options['max_iter']=3000
            # ipopt.options['acceptable_tol']=1e-6
            unit.initialize()
            unit.display()
            return
        elif isinstance(unit,idaes_pressure_changer.PressureChanger):
            unit.display()
            unit.initialize()
            unit.costing.initialize()
            ipopt.options['max_iter']=300
            ipopt.options['acceptable_tol']=1e-2
            scaled_unit = pyo.TransformationFactory('core.scale_model').create_using(unit)
            ipopt.solve(scaled_unit,tee=False)
            pyo.TransformationFactory('core.scale_model').propagate_solution(scaled_unit,unit)
            ipopt.options['max_iter']=3000
            ipopt.options['acceptable_tol']=1e-6
            unit.display()
            return
        else:
            unit.initialize()
            return

    seq.run(m,custom_initializer_func)

    #undo initialization dof
    m.fs.source_comp.outlet.pressure[0].unfix()
    m.fs.recycle_comp.outlet.pressure[0].unfix()
    m.fs.purge_splitter.split_fraction[0,'to_purge'].unfix()
    m.fs.purge_splitter.to_purge.pressure[0].unfix()
    #fix pressure drops instead of diameters
    for i in range(0,15):
        pipe = getattr(m.fs,"pipe"+str(i))
        pipe.Pdrop_constraint.deactivate()
        pipe.max_pressure.unfix()

    #splitter outlet pressures
    for j in m.fs.split0.create_outlet_list():
        getattr(m.fs.split0,j).pressure[0].unfix()
    for j in m.fs.split1.create_outlet_list():
        getattr(m.fs.split1,j).pressure[0].unfix()
    for j in m.fs.split2.create_outlet_list():
        getattr(m.fs.split2,j).pressure[0].unfix()
    for j in m.fs.split3.create_outlet_list():
        getattr(m.fs.split3,j).pressure[0].unfix()

    #mixer outlet pressures
    m.fs.mix0.outlet.pressure[0].unfix()
    m.fs.mix1.outlet.pressure[0].unfix()
    m.fs.mix2.outlet.pressure[0].unfix()
    m.fs.mix3.outlet.pressure[0].unfix()
    m.fs.mixPF.outlet.pressure[0].unfix()

    #other compressors
    for i in range(1,7):
        getattr(m.fs, "comp"+str(i)).outlet.pressure[0].unfix()

    #processing facility d.o.f
    m.fs.processing_facility.loop_solvent_flow.unfix()
    m.fs.processing_facility.exchanger1.outlet_temperature.unfix()
    m.fs.processing_facility.flash3.pressure.unfix()
    m.fs.processing_facility.column.Nstages.unfix()
    m.fs.processing_facility.activate_surrogate_input_bounds()

    #guess inlet flowrate
    m.fs.source_comp.inlet.flow_mol[0,'co2'].unfix()

from co2_eor.utility.util_funcs import run_with_redirection
run_with_redirection(initialization1,output_dir+m_name+'_initialization1_log.txt')
print(f'initialization 1 done')
assert idaes.core.util.model_statistics.degrees_of_freedom(m) == DoF

#build economics
m.fs.co2_price = pyo.Param(initialize=0.068) #$/kg
m.fs.electricity_price = pyo.Param(initialize=5E-8) #$/J
m.fs.refrigeration_price = pyo.Param(initialize=3.9E-8) #$/J
m.fs.capex_recovery_factor = pyo.Param(initialize=0.1315)

@m.fs.Expression()
def capex(fs):
    expr = fs.source_comp.costing.capital_cost + fs.recycle_comp.costing.capital_cost + fs.processing_facility.costing.capital_cost
    for i in range(1,7):
        expr += getattr(fs,"comp"+str(i)).costing.capital_cost
    for i in range(0,15):
        expr += getattr(fs,"pipe"+str(i)).costing.capital_cost
    return expr

@m.fs.Expression()
def opex(fs):
    electrical_power = fs.source_comp.work_mechanical[0] + fs.recycle_comp.work_mechanical[0] + m.fs.processing_facility.costing.electrical_power
    for i in range(1,7):
        electrical_power += getattr(fs,"comp"+str(i)).work_mechanical[0]
    return 365*24*3600* (m.fs.electricity_price*electrical_power+ m.fs.processing_facility.costing.refrigeration_duty*m.fs.refrigeration_price)

m.fs.raw_mats = pyo.Expression(expr=m.fs.source_comp.control_volume.properties_in[0].flow_mass*m.fs.co2_price*3600*24*365)

@m.fs.Expression()
def total_oil_production(fs):
    expr = 0
    for i in range(1,7):
        expr += getattr(fs,"well"+str(i)).q_OIL_PROD
    return expr
m.fs.production_target = pyo.Constraint(expr=m.fs.total_oil_production *3600*24*365 >= 450000)
set_scaling_factor(m.fs.production_target,1e-6)

m.fs.obj = pyo.Objective(expr=
                         m.fs.capex_recovery_factor*m.fs.capex
                         +m.fs.opex
                         +m.fs.raw_mats)
set_scaling_factor(m.fs.obj,1e-7)

with open(output_dir+m_name+'_postinitialization1_display.txt','w') as f:
    with contextlib.redirect_stdout(f):
        m.display()

from co2_eor.utility.util_funcs import report_infeasibilities
report_infeasibilities(m, output_dir+m_name+"_presolve_violations.log")

from idaes.core.scaling.autoscaling import AutoScaler
autoScaler=AutoScaler(overwrite=True)
autoScaler.scale_variables_by_magnitude(m)
#autoScaler2=AutoScaler(overwrite=True)
autoScaler.scale_constraints_by_jacobian_norm(m)

m.write(output_dir+m_name+"_model.nl", format="nl")

ipopt.options['max_iter'] = 30000
ipopt.options['linear_solver']='ma27'
ipopt.options['acceptable_tol']=1E-6
ipopt.options['tol']=1E-8
#ipopt.options['nlp_scaling_method']='user-scaling'
ipopt.options['expect_infeasible_problem']='yes'
#ipopt.options['mu_init'] = 1
#ipopt.options['OF_evaluate_orig_obj_at_resto_trial']='no'
ipopt.options['OF_start_with_resto']='yes'
#ipopt.options['OF_soft_resto_pderror_reduction_factor']=0
#ipopt.options['OF_required_infeasibility_reduction']=0.99
#scale model
scaled_m = pyo.TransformationFactory("core.scale_model").create_using(m)
res=ipopt.solve(scaled_m,tee=True,keepfiles=True,symbolic_solver_labels=True,)
pyo.TransformationFactory("core.scale_model").propagate_solution(scaled_m,m)

with open(m.fs, output_dir+m_name+'_postsolve_display.txt', 'w') as f:
    with contextlib.redirect_stdout(f):
        m.display()

report_infeasibilities(m, output_dir+m_name+"_postsolve_violations.log")

pyo.assert_optimal_termination(res)

from co2_eor.utility.util_funcs import export_flowsheet_to_excel
export_flowsheet_to_excel(m.fs, output_dir+m_name+'_solve.xlsx')

m.fs.visualize(m_name)
time.sleep(5)
