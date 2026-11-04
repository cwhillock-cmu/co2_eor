import pyomo.environ as pyo
import idaes.core
from pyomo.network import Arc
from pyomo.environ import units
import contextlib
import numpy as np

import idaes.models.properties.general_helmholtz as idaesHelmholtz
import idaes.models.unit_models.pressure_changer as idaesPressureChanger
from idaes.core.util.scaling import set_scaling_factor

from idaes.core import MomentumBalanceType
from co2_eor.splitter_unit import EnergySplittingType
from co2_eor.mixer_unit import MomentumMixingType

from co2_eor.SSLW import SSLWCosting, SSLWCostingData
from co2_eor.SSLW import CompressorType, CompressorDriveType, CompressorMaterial
from co2_eor.SSLW import VesselMaterial

from co2_eor import mixer,splitter, heater

from idaes.models.properties.modular_properties.base.generic_property import GenericParameterBlock
from co2_eor.liq_pipe import liqPipe
from co2_eor.wellpattern_unit import wellpattern
import co2_eor.MPF.thermo_config as thermo_config
from co2_eor.gas_pipe import gasPipe
from co2_eor.processing_facility import processingFacility
from co2_eor.surrogates.rectisol_surrogate import rectisolSurrogate
from co2_eor.mpf_to_helmholtz import mpf_helmholtz_converter

from co2_eor.utility.solvers import ipopt, conopt

filename = "bakken_network_rev3_surrogate_cubic"

#model, flowsheet, and parameter block
m = pyo.ConcreteModel()
m.fs = idaes.core.FlowsheetBlock(dynamic=False)
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

compressor_config = {
        "property_package": m.fs.props_mix_liq_cubic,
        "dynamic":False,
        "compressor":True,
        "thermodynamic_assumption":idaesPressureChanger.ThermodynamicAssumption.isentropic,
        }

liqPipe_config = {
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
        "momentum_mixing_type":MomentumMixingType.inequality
        }

splitter_config = {
        "property_package":m.fs.props_mix_liq_cubic,
        "ideal_separation":False,
        "energy_split_basis":EnergySplittingType.equal_molar_enthalpy,
        "momentum_balance_type":MomentumBalanceType.none,
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
    "injectivity_index":2.93e-10,#2.93e-10,
    "reference_pressure":413*100000,
    "use_correction_factor":False,
    "depth":1000,
    "pressure_max":650*100000, #650?
    "outlet_pressure":30*100000,
    "outlet_temperature":300,
    "multiplier":1
}

heater_config = {
    "property_package":m.fs.props_mix_liq_cubic,
}

gasPipe_config = {
    "property_package":m.fs.props_mix_vap_cubic,
    "average_pressure_type":'nonlinear',
    "allowable_stress":1380*100000,
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

pipe_costing_config = {
    "flowsheet_costing_block":m.fs.costing,
    "costing_method":SSLWCostingData.cost_horizontal_vessel,
    "costing_method_arguments":{
        "material_type":VesselMaterial.StainlessSteel316,
        "include_platforms_ladders":False,
    }
}

#prob not needed
gasPipe_costing_config = {
    "flowsheet_costing_block":m.fs.costing,
    "costing_method":SSLWCostingData.cost_horizontal_vessel,
    "costing_method_arguments":{
        "material_type":VesselMaterial.StainlessSteel316,
        "include_platforms_ladders":False,
    }
}

surrogate_input_bounds = [
    [400,200,50e5,320,2000,260,6e5,10.5],
    [100,50,10e5,290,800,220,2e5,3.5]
]
hyperplane_list = [
    {"weights": [1.12241413e-03, -2.56485953e-02, 1.67543279e-06, 2.40731609e-03, 1.42386045e-03, 1.04199595e-02, -6.81482689e-06, 9.81953962e-01],
     "bias":-15.0399895}
]
surrogate_object = idaes.core.surrogate.AlamoSurrogate.load_from_file('temps/alamo_surrogate.json')
processingFacility_config = {
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
        "recycle_outlet_temperature":350
    },
    "surrogate_object":surrogate_object,
    "hyperplane_cuts":hyperplane_list,
}

#distribution network:
#source branch
# m.fs.source_comp = idaesPressureChanger.PressureChanger(**compressor_config)
# m.fs.source_comp.costing = idaes.core.UnitModelCostingBlock(**compressor_costing_config)
# m.fs.source_comp.efficiency_isentropic.fix(0.85)

m.fs.pipe0 = liqPipe(**liqPipe_config,length=15000)
m.fs.pipe0.costing = idaes.core.UnitModelCostingBlock(**pipe_costing_config)
m.fs.pipe0.diameter = 1.2
m.fs.pipe0.roughness.fix(0.0475e-3)
m.fs.pipe0.max_pressure.fix(650*100000)

m.fs.recycle_comp = idaesPressureChanger.PressureChanger(**compressor_config)
m.fs.recycle_comp.costing = idaes.core.UnitModelCostingBlock(**compressor_costing_config)
m.fs.recycle_comp.efficiency_isentropic.fix(0.85)
m.fs.recycle_comp.pressure_minimum_eq = pyo.Constraint(expr=m.fs.recycle_comp.outlet.pressure[0]>=80*100000)

m.fs.chiller = heater(**heater_config)
#m.fs.chiller.outlet_temp_eq = pyo.Constraint(expr=m.fs.chiller.control_volume.properties_out[0].temperature<=320)
#m.fs.chiller.outlet_temp_eq2 = pyo.Constraint(expr=m.fs.chiller.control_volume.properties_out[0].temperature>=290)
#m.fs.chiller.positive_flow_eq = pyo.Constraint(expr=m.fs.chiller.control_volume.properties_in[0].flow_mol>=0)
#m.fs.chiller.costing = idaes.core.UnitModelCostingBlock(
#    flowsheet_costing_block = m.fs.costing,
#    costing_method = SSLWCostingData.cost_heat_exchanger,
#)
#m.fs.chiller.control_volume.properties_out[0].temperature.fix(320)
m.fs.chiller.outlet_temp_eq = pyo.Constraint(expr=m.fs.chiller.control_volume.properties_out[0].temperature==318)

#m.fs.s_source_comp_pipe0 = Arc(source=m.fs.source_comp.outlet,destination=m.fs.pipe0.inlet)

m.fs.mix0 = mixer(**mixer_config,property_package=m.fs.props_mix_liq_cubic,inlet_list=['from_pipe0','from_recycle_comp'])

#m.fs.s_recycle_comp_chiller = Arc(source=m.fs.recycle_comp.outlet,destination=m.fs.chiller.inlet)
#m.fs.s_chiller_mix0 = Arc(source=m.fs.chiller.outlet,destination=m.fs.mix0.from_chiller)
m.fs.s_recycle_comp_mix0 = Arc(source=m.fs.recycle_comp.outlet,destination=m.fs.mix0.from_recycle_comp)
m.fs.s_pipe0_mix0 = Arc(source=m.fs.pipe0.outlet,destination=m.fs.mix0.from_pipe0)
m.fs.s_mix0_chiller = Arc(source=m.fs.mix0.outlet,destination=m.fs.chiller.inlet)

m.fs.split0 = splitter(**splitter_config,outlet_list=['to_pipe1','to_pipe3','to_pipe6'])
#m.fs.s_mix0_split0 = Arc(source=m.fs.mix0.outlet,destination=m.fs.split0.inlet)
m.fs.s_chiller_split0 = Arc(source=m.fs.chiller.outlet,destination=m.fs.split0.inlet)

#first branch, wells 1 and 2
m.fs.pipe1 = liqPipe(**liqPipe_config,length=800)
m.fs.pipe1.costing=idaes.core.UnitModelCostingBlock(**pipe_costing_config)
#m.fs.pipe1.diameter.fix(0.5)
m.fs.pipe1.roughness.fix(0.0475e-3)
m.fs.pipe1.max_pressure.fix(650*100000)

m.fs.split1 = splitter(**splitter_config,outlet_list=['to_comp1','to_pipe2'])

m.fs.comp1 = idaesPressureChanger.PressureChanger(**compressor_config)
m.fs.comp1.costing = idaes.core.UnitModelCostingBlock(**compressor_costing_config)
m.fs.comp1.efficiency_isentropic.fix(0.85)
m.fs.comp1.pressure_minimum_eq = pyo.Constraint(expr=m.fs.comp1.inlet.pressure[0]>=80*100000)

m.fs.well1 = wellpattern(**wellpattern_config)
m.fs.well1.HCPV.fix(1.1)

m.fs.pipe2 = liqPipe(**liqPipe_config,length=1600)
m.fs.pipe2.costing=idaes.core.UnitModelCostingBlock(**pipe_costing_config)
#m.fs.pipe2.diameter.fix(0.5)
m.fs.pipe2.roughness.fix(0.0475e-3)
m.fs.pipe2.max_pressure.fix(650*100000)

m.fs.comp2 = idaesPressureChanger.PressureChanger(**compressor_config)
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
m.fs.pipe3 = liqPipe(**liqPipe_config,length=3200)
m.fs.pipe3.costing=idaes.core.UnitModelCostingBlock(**pipe_costing_config)
#m.fs.pipe3.diameter.fix(0.5)
m.fs.pipe3.roughness.fix(0.0475e-3)
m.fs.pipe3.max_pressure.fix(650*100000)

m.fs.split2 = splitter(**splitter_config,outlet_list=['to_pipe4','to_pipe5'])

m.fs.pipe4 = liqPipe(**liqPipe_config,length=4000)
m.fs.pipe4.costing=idaes.core.UnitModelCostingBlock(**pipe_costing_config)
#m.fs.pipe4.diameter.fix(0.5)
m.fs.pipe4.roughness.fix(0.0475e-3)
m.fs.pipe4.max_pressure.fix(650*100000)

m.fs.comp3 = idaesPressureChanger.PressureChanger(**compressor_config)
m.fs.comp3.costing = idaes.core.UnitModelCostingBlock(**compressor_costing_config)
m.fs.comp3.efficiency_isentropic.fix(0.85)
m.fs.comp3.pressure_minimum_eq = pyo.Constraint(expr=m.fs.comp3.inlet.pressure[0]>=80*100000)

m.fs.well3 = wellpattern(**wellpattern_config)
m.fs.well3.HCPV.fix(0.8)

m.fs.pipe5 = liqPipe(**liqPipe_config,length=2400)
m.fs.pipe5.costing=idaes.core.UnitModelCostingBlock(**pipe_costing_config)
#m.fs.pipe5.diameter.fix(1)
m.fs.pipe5.roughness.fix(0.0475e-3)
m.fs.pipe5.max_pressure.fix(650*100000)

m.fs.comp4 = idaesPressureChanger.PressureChanger(**compressor_config)
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
m.fs.pipe6 = liqPipe(**liqPipe_config,length=3200)
m.fs.pipe6.costing=idaes.core.UnitModelCostingBlock(**pipe_costing_config)
#m.fs.pipe6.diameter.fix(0.5)
m.fs.pipe6.roughness.fix(0.0475e-3)
m.fs.pipe6.max_pressure.fix(650*100000)

m.fs.split3 = splitter(**splitter_config,outlet_list=['to_pipe7','to_pipe8'])

#m.fs.pipe7 = liqPipe(**liqPipe_config,length=2400)
m.fs.pipe7 = liqPipe(**liqPipe_config,length=2400)
m.fs.pipe7.costing=idaes.core.UnitModelCostingBlock(**pipe_costing_config)
#m.fs.pipe7.diameter.fix(0.9)
m.fs.pipe7.roughness.fix(0.0475e-3)
m.fs.pipe7.max_pressure.fix(650*100000)

m.fs.comp5 = idaesPressureChanger.PressureChanger(**compressor_config)
m.fs.comp5.costing = idaes.core.UnitModelCostingBlock(**compressor_costing_config)
m.fs.comp5.efficiency_isentropic.fix(0.85)
m.fs.comp5.pressure_minimum_eq = pyo.Constraint(expr=m.fs.comp5.inlet.pressure[0]>=80*100000)

m.fs.well5 = wellpattern(**wellpattern_config)
m.fs.well5.HCPV.fix(1.0)

m.fs.pipe8 = liqPipe(**liqPipe_config,length=4000)
m.fs.pipe8.costing=idaes.core.UnitModelCostingBlock(**pipe_costing_config)
#m.fs.pipe8.diameter.fix(1)
m.fs.pipe8.roughness.fix(0.0475e-3)
m.fs.pipe8.max_pressure.fix(650*100000)

m.fs.comp6 = idaesPressureChanger.PressureChanger(**compressor_config)
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
m.fs.pipe9 = gasPipe(**gasPipe_config,length=1000)
m.fs.pipe9.costing=idaes.core.UnitModelCostingBlock(**pipe_costing_config)
m.fs.pipe9.diameter = 0.6
m.fs.pipe9.roughness.fix(0.0475e-3)
m.fs.pipe9.max_pressure.fix(50*100000)

m.fs.mix1 = mixer(**mixer_config,property_package=m.fs.props_mix_vap_cubic,inlet_list=['from_well1','from_pipe10'])

m.fs.pipe10 = gasPipe(**gasPipe_config,length=1600)
m.fs.pipe10.costing=idaes.core.UnitModelCostingBlock(**pipe_costing_config)
m.fs.pipe10.diameter = 0.6
m.fs.pipe10.roughness.fix(0.0475e-3)
m.fs.pipe10.max_pressure.fix(50*100000)

m.fs.mix2 = mixer(**mixer_config,property_package=m.fs.props_mix_vap_cubic,inlet_list=['from_well2','from_pipe11'])

m.fs.pipe11 = gasPipe(**gasPipe_config,length=3200)
m.fs.pipe11.costing=idaes.core.UnitModelCostingBlock(**pipe_costing_config)
m.fs.pipe11.diameter = 0.6
m.fs.pipe11.roughness.fix(0.0475e-3)
m.fs.pipe11.max_pressure.fix(50*100000)

m.fs.s_well4_pipe11 = Arc(source=m.fs.well4.outlet,destination=m.fs.pipe11.inlet)
m.fs.s_pipe11_mix2 = Arc(source=m.fs.pipe11.outlet,destination=m.fs.mix2.from_pipe11)
m.fs.s_well2_mix2 = Arc(source=m.fs.well2.outlet,destination=m.fs.mix2.from_well2)
m.fs.s_mix2_pipe10 = Arc(source=m.fs.mix2.outlet,destination=m.fs.pipe10.inlet)
m.fs.s_pipe10_mix1 = Arc(source=m.fs.pipe10.outlet,destination=m.fs.mix1.from_pipe10)
m.fs.s_well1_mix1 = Arc(source=m.fs.well1.outlet,destination=m.fs.mix1.from_well1)
m.fs.s_mix1_pipe9 = Arc(source=m.fs.mix1.outlet,destination=m.fs.pipe9.inlet)
m.fs.s_pipe9_mixPF = Arc(source=m.fs.pipe9.outlet,destination=m.fs.mixPF.from_pipe9)

#collection from wells 5,6
m.fs.pipe12 = gasPipe(**gasPipe_config,length=2000)
m.fs.pipe12.costing=idaes.core.UnitModelCostingBlock(**pipe_costing_config)
m.fs.pipe12.diameter = 0.6
m.fs.pipe12.roughness.fix(0.0475e-3)
m.fs.pipe12.max_pressure.fix(50*100000)

m.fs.mix3 = mixer(**mixer_config,property_package=m.fs.props_mix_vap_cubic,inlet_list=['from_well6','from_pipe13'])

#m.fs.pipe13 = gasPipe(**gasPipe_config,length=5500)
m.fs.pipe13 = gasPipe(**gasPipe_config,length=5500)
m.fs.pipe13.costing=idaes.core.UnitModelCostingBlock(**pipe_costing_config)
m.fs.pipe13.diameter = 0.6
m.fs.pipe13.roughness.fix(0.0475e-3)
m.fs.pipe13.max_pressure.fix(50*100000)

m.fs.s_well5_pipe13 = Arc(source=m.fs.well5.outlet,destination=m.fs.pipe13.inlet)
m.fs.s_pipe13_mix3 = Arc(source=m.fs.pipe13.outlet,destination=m.fs.mix3.from_pipe13)
m.fs.s_well6_mix3 = Arc(source=m.fs.well6.outlet,destination=m.fs.mix3.from_well6)
m.fs.s_mix3_pipe12 = Arc(source=m.fs.mix3.outlet,destination=m.fs.pipe12.inlet)
m.fs.s_pipe12_mixPF = Arc(source=m.fs.pipe12.outlet,destination=m.fs.mixPF.from_pipe12)

#collection from well 3
m.fs.pipe14 = gasPipe(**gasPipe_config,length=2800)
m.fs.pipe14.costing=idaes.core.UnitModelCostingBlock(**pipe_costing_config)
m.fs.pipe14.diameter = 0.6
m.fs.pipe14.roughness.fix(0.0475e-3)
m.fs.pipe14.max_pressure.fix(50*100000)

m.fs.s_well3_pipe14 = Arc(source=m.fs.well3.outlet,destination=m.fs.pipe14.inlet)
m.fs.s_pipe14_mixPF = Arc(source=m.fs.pipe14.outlet,destination=m.fs.mixPF.from_pipe14)

#processing facility to recycle compressor
#m.fs.processing_facility=rectisolSurrogate(**processingFacility_config)
m.fs.processing_facility = processingFacility(
    inlet_property_package=m.fs.props_mix_vap_cubic,
    outlet_property_package=m.fs.props_mix_vap_cubic,
    recycle_property_package=m.fs.props_mix_liq_cubic,
    key_component='co2',
    minimum_inlet_pressure=10*100000,
    recycle_pressure=74*100000,
    recycle_temperature=350,
    outlet_pressure=5*100000
)

m.fs.purge_splitter = splitter(**splitter_config,outlet_list=['purge','to_recycle_comp'])
#m.fs.purge_splitter.positive_purge_eq = pyo.Constraint(expr=m.fs.purge_splitter.purge.flow_mol[0]>=0) #enforce material does not go into system from purge stream
m.fs.purge_splitter.purge.flow_mol_comp[0,'co2'].lb = 0
m.fs.purge_splitter.purge.flow_mol_comp[0,'ch4'].lb = 0

m.fs.s_mixPF_processing_facility = Arc(source=m.fs.mixPF.outlet,destination=m.fs.processing_facility.inlet)
#m.fs.s_procFac_eos_converter = Arc(source=m.fs.procFac.recycle,destination=m.fs.mpfHelmholtzConverter.inlet)
#m.fs.s_eos_converter_purge_splitter = Arc(source=m.fs.mpfHelmholtzConverter.outlet,destination=m.fs.purge_splitter.inlet)
m.fs.s_processing_facility_purge_splitter = Arc(source=m.fs.processing_facility.recycle,destination=m.fs.purge_splitter.inlet)
m.fs.s_purge_splitter_recycle_comp = Arc(source=m.fs.purge_splitter.to_recycle_comp,destination=m.fs.recycle_comp.inlet)

#expand arcs
pyo.TransformationFactory("network.expand_arcs").apply_to(m)

#inlet specifications
# m.fs.source_comp.inlet.pressure[0].fix(100*100000)
# #m.fs.source_comp.inlet.enth_mol[0].fix(m.fs.props1.htpx(T=300*units.K,p=100*100000*units.Pa,amount_basis=idaesHelmholtz.AmountBasis.MOLE))
# m.fs.source_comp.inlet.temperature[0].fix(300)
m.fs.pipe0.inlet.pressure[0].fix(201*100000)
m.fs.pipe0.inlet.temperature[0].fix(300)

with open('temps/'+filename+'_preinitialization1_pprint.txt', 'w') as f:
    with contextlib.redirect_stdout(f):
        m.pprint()

#fix initialization 1 dof

for i in range(0,9):
    pipe = getattr(m.fs,"pipe"+str(i))
    pipe.Pdrop_constraint = pyo.Constraint(expr=pipe.Pdrop==1*100000)
    set_scaling_factor(pipe.Pdrop_constraint,1e-5)
for i in range(9,15):
    pipe = getattr(m.fs,"pipe"+str(i))
    pipe.Pdrop_constraint = pyo.Constraint(expr=pipe.Pdrop==5*100000)
    set_scaling_factor(pipe.Pdrop_constraint,1e-5)
m.fs.split0.to_pipe1.pressure[0].fix(195*100000)
m.fs.split0.to_pipe3.pressure[0].fix(195*100000)
m.fs.split0.to_pipe6.pressure[0].fix(195*100000)
m.fs.split1.to_comp1.pressure[0].fix(190*100000)
m.fs.split1.to_pipe2.pressure[0].fix(190*100000)
m.fs.split2.to_pipe4.pressure[0].fix(185*100000)
m.fs.split2.to_pipe5.pressure[0].fix(185*100000)
m.fs.split3.to_pipe7.pressure[0].fix(185*100000)
m.fs.split3.to_pipe8.pressure[0].fix(185*100000)

m.fs.mix1.outlet.pressure[0].fix(15*100000)
m.fs.mix2.outlet.pressure[0].fix(25*100000)
m.fs.mix3.outlet.pressure[0].fix(25*100000)
m.fs.mixPF.outlet.pressure[0].fix(10*100000)

# m.fs.source_comp.outlet.pressure[0].fix(200*100000)
# m.fs.source_comp.inlet.flow_mol_comp[0,'ch4'].fix(1e-6)
m.fs.pipe0.inlet.flow_mol_comp[0,'ch4'].fix(30)

m.fs.recycle_comp.outlet.pressure[0].fix(200*100000)
m.fs.mix0.outlet.pressure[0].fix(200*100000)

for i in range(1,7):
    getattr(m.fs, "comp"+str(i)).outlet.pressure[0].fix(550*100000)

m.fs.purge_splitter.split_fraction[0,'purge'].fix(1e-14)
#m.fs.purge_splitter.purge.pressure[0].fix(73.9*100000)
m.fs.purge_splitter.to_recycle_comp.pressure[0].fix(73.9*100000)

#m.fs.source_comp.inlet.flow_mol_comp[0,'co2'].fix(160)
m.fs.pipe0.inlet.flow_mol_comp[0,'co2'].fix(180)

#check for square model (outside of processing facility d.o.f)
print(f'number of variables={len(list(m.component_data_objects(pyo.Var)))}')
print(f'number of constraints={len(list(m.component_data_objects(pyo.Constraint)))}')
DoF = idaes.core.util.model_statistics.degrees_of_freedom(m)
print(f'degrees of freedom={DoF}')
#input('paused')

print(f'starting solve')

def initialization_type1(m):
    from pyomo.network import SequentialDecomposition
    from idaes.core.util.initialization import propagate_state
    seq = SequentialDecomposition(tol=1e-6,tol_type='rel')
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
    # tear_guesses2 = {
    #     "flow_mol":{
    #         (0,'co2'):250,
    #         (0,'ch4'):8,
    #     },
    #     "temperature":{0:449},
    #     #"enth_mol":{0:m.fs.props1.htpx(T=300*units.K,p=200*100000*units.Pa,amount_basis=idaesHelmholtz.AmountBasis.MOLE)},
    #     "pressure":{0:200*100000},
    # }
    tear_guesses2 = {
        "flow_mol":{
            (0,'co2'):355,
            (0,'ch4'):30.2,
        },
        "temperature":{0:448},
        #"enth_mol":{0:m.fs.props1.htpx(T=300*units.K,p=200*100000*units.Pa,amount_basis=idaesHelmholtz.AmountBasis.MOLE)},
        "pressure":{0:200*100000},
    }
    #seq.set_guesses_for(m.fs.recycle_comp.inlet,tear_guesses1)
    seq.set_guesses_for(m.fs.mix0.from_recycle_comp,tear_guesses2)

    def custom_initializer_new(unit):
        #unit.display()
        if isinstance(unit,(liqPipe,gasPipe)):
            #unit.display()
            unit.initialize(tee=False,display_after=True)
            #input('paused')
            return
        elif isinstance(unit,wellpattern):
            #unit.display()
            #unit.initialize(tee=False,display_after=True)
            #input('paused')
            number = unit.name[-1]
            comp = getattr(m.fs,"comp"+number)
            unit.inlet.pressure[0].unfix()
            unit.initialize()
            comp.outlet.pressure[0].fix(pyo.value(unit.inlet.pressure[0]))
            comp.initialize()
            comp.costing.initialize()
            return
        elif unit.name == 'fs.processing_facility':
            #unit.display()
            unit.initialize(display_after=True)
            #unit.custom_propagate_state()
            #unit.display()
            #input('paused')
            return
        elif isinstance(unit,mixer):
            unit.display()
            # ipopt.options['max_iter']=30
            # ipopt.options['acceptable_tol']=1e-2
            set_scaling_factor(unit.enthalpy_mixing_equations,1e-2)
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
            unit.display()
            ipopt.options['max_iter']=300
            ipopt.options['acceptable_tol']=1e-2
            set_scaling_factor(unit.molar_enthalpy_equality_eqn,1e-3)
            scaled_unit = pyo.TransformationFactory('core.scale_model').create_using(unit)
            ipopt.solve(scaled_unit,tee=False)
            pyo.TransformationFactory('core.scale_model').propagate_solution(scaled_unit,unit)
            ipopt.options['max_iter']=3000
            ipopt.options['acceptable_tol']=1e-6
            unit.initialize()
            #unit.display()
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
for i in range(1,7):
    unit = getattr(m.fs,"comp"+str(i))
    unit.inlet.flow_mol_comp[0,'co2'].fix()
    unit.inlet.flow_mol_comp[0,'ch4'].fix()
    unit.inlet.pressure[0].fix()
    unit.inlet.temperature[0].fix()
    ipopt.options['max_iter']=300
    ipopt.options['acceptable_tol']=1e-2
    scaled_unit = pyo.TransformationFactory('core.scale_model').create_using(unit)
    ipopt.solve(scaled_unit,tee=False)
    pyo.TransformationFactory('core.scale_model').propagate_solution(scaled_unit,unit)
    ipopt.options['max_iter']=3000
    ipopt.options['acceptable_tol']=1e-6
    unit.inlet.flow_mol_comp[0,'co2'].unfix()
    unit.inlet.flow_mol_comp[0,'ch4'].unfix()
    unit.inlet.pressure[0].unfix()
    unit.inlet.temperature[0].unfix()
print(f'initialization 1 done')

# #unfix initialization 1 degrees of freedom
# #m.fs.source_comp.inlet.flow_mol_comp[0,'co2'].unfix()
m.fs.pipe0.inlet.flow_mol_comp[0,'co2'].unfix()
# #m.fs.source_comp.inlet.flow_mol_comp[0,'ch4'].fix(0)
# #m.fs.source_comp.outlet.pressure[0].unfix()
# m.fs.recycle_comp.outlet.pressure[0].unfix()
# m.fs.mix0.outlet.pressure[0].unfix()
# m.fs.purge_splitter.to_recycle_comp.pressure[0].unfix()

# for i in range(1,7):
#     getattr(m.fs, "comp"+str(i)).outlet.pressure[0].unfix()

# for i in range(0,15):
#     pipe = getattr(m.fs,"pipe"+str(i))
#     pipe.Pdrop_constraint.deactivate()
# m.fs.split0.to_pipe1.pressure[0].unfix()
# m.fs.split0.to_pipe3.pressure[0].unfix()
# m.fs.split0.to_pipe6.pressure[0].unfix()
# m.fs.split1.to_comp1.pressure[0].unfix()
# m.fs.split1.to_pipe2.pressure[0].unfix()
# m.fs.split2.to_pipe4.pressure[0].unfix()
# m.fs.split2.to_pipe5.pressure[0].unfix()
# m.fs.split3.to_pipe7.pressure[0].unfix()
# m.fs.split3.to_pipe8.pressure[0].unfix()
# #m.fs.purge_splitter.split_fraction[0,'purge'].unfix()

# m.fs.mix1.outlet.pressure[0].unfix()
# m.fs.mix2.outlet.pressure[0].unfix()
# m.fs.mix3.outlet.pressure[0].unfix()
# m.fs.mixPF.outlet.pressure[0].unfix()

#create production target
m.fs.prod_target = pyo.Param(mutable=True,initialize=400000) #STB/yr
#production target constraint
@m.fs.Constraint()
def production_target_constraint(fs):
    rhs = 0
    for i in range(1,7):
        rhs += getattr(fs,"well"+str(i)).q_OIL_PROD
    return fs.prod_target / 365 <= rhs * 3600 * 24
#m.fs.production_target_constraint.deactivate()

m.fs.CW_price = pyo.Param(initialize=0.06) #$/m3
m.fs.electricity_price = pyo.Param(initialize=5E-8) #$/J
m.fs.refrigeration_price = pyo.Param(initialize=3.9E-8) #$/J
#create economic terms
@m.fs.Expression()
def opex(fs):
    electrical_power = (0+#fs.source_comp.work_mechanical[0] 
                        + fs.recycle_comp.work_mechanical[0] 
                        #+ fs.processing_facility.costing.electrical_power
                        # #- 0.2*fs.chiller.heat_duty[0]
    )
    for i in range(1,7):
        electrical_power += getattr(fs,"comp"+str(i)).work_mechanical[0]
    return (365*24*3600* (fs.electricity_price*electrical_power 
                          #+ fs.refrigeration_price*fs.processing_facility.costing.refrigeration_duty
                          )
                          + fs.processing_facility.costing.opex
                          )

m.fs.co2_price = pyo.Param(initialize=0.068)
m.fs.raw_mats = pyo.Expression(
    expr=365*24*3600* m.fs.co2_price*m.fs.pipe0.control_volume.properties_in[0].flow_mass
)

m.fs.oil_price = pyo.Param(initialize=700,mutable=True)
@m.fs.Expression()
def revenue(fs):
    oil_prod = 0
    for i in range(1,7):
        oil_prod += getattr(fs,"well"+str(i)).q_OIL_PROD
    return 365*24*3600* fs.oil_price*oil_prod

@m.fs.Expression()
def capex(fs):
    expr = (0#+fs.source_comp.costing.capital_cost 
            + fs.recycle_comp.costing.capital_cost 
            + fs.processing_facility.costing.capital_cost
            # + fs.chiller.costing.capital_cost
    )
    for i in range(1,7):
        expr += getattr(fs,"comp"+str(i)).costing.capital_cost
    for i in range(0,15):
        expr += getattr(fs,"pipe"+str(i)).costing.capital_cost
    return expr

m.fs.capital_recovery_factor = pyo.Param(initialize=0.1315)
m.fs.obj_no_revenue = pyo.Objective(
    expr=(m.fs.capex*m.fs.capital_recovery_factor+m.fs.opex+m.fs.raw_mats)
)
set_scaling_factor(m.fs.obj_no_revenue,1e-7)
m.fs.test_obj_constraint = pyo.Constraint(expr=(m.fs.capex*m.fs.capital_recovery_factor+m.fs.opex+m.fs.raw_mats)>=0)

with open('temps/'+filename+'_presolve_display.txt', 'w') as f:
    with contextlib.redirect_stdout(f):
        m.display()

from co2_eor.utility.util_funcs import report_infeasibilities
report_infeasibilities(m, "temps/presolve_violations.log")

from idaes.core.scaling.autoscaling import AutoScaler
autoScaler=AutoScaler(overwrite=False)
autoScaler.scale_variables_by_magnitude(m)
autoScaler.scale_constraints_by_jacobian_norm(m)

#scale model
scaled_m = pyo.TransformationFactory("core.scale_model").create_using(m)
# with open('temps/'+filename+'_scaled_presolve_display.txt', 'w') as f:
#     with contextlib.redirect_stdout(f):
#         scaled_m.display()
report_infeasibilities(scaled_m, "temps/scaled_presolve_violations.log")
#solve flowsheet
ipopt.options['max_iter'] = 2000
ipopt.options['linear_solver']='ma27'
ipopt.options['acceptable_tol']=1E-6
ipopt.options['tol']=1E-8
ipopt.options['nlp_scaling_method']='gradient-based'
#ipopt.options['nlp_scaling_method']='user-scaling'
ipopt.options['expect_infeasible_problem']='yes'
#ipopt.options['OF_evaluate_orig_obj_at_resto_trial']='no'
# ipopt.options['OF_start_with_resto']='yes'
# ipopt.options['OF_soft_resto_pderror_reduction_factor']=0
# ipopt.options['OF_required_infeasibility_reduction']=0.01
# ipopt.options['OF_obj_max_inc']=10
# scaled_m.fs.scaled_obj_no_revenue.deactivate()
# res=ipopt.solve(scaled_m,tee=True)
# ipopt.options['linear_solver']='ma97'
# ipopt.options['expect_infeasible_problem']='no'
# scaled_m.fs.scaled_obj_no_revenue.activate()
res=ipopt.solve(scaled_m,tee=True)
report_infeasibilities(scaled_m, "temps/scaled_postsolve_violations.log")
#unscale model
pyo.TransformationFactory("core.scale_model").propagate_solution(scaled_m,m)

with open('temps/'+filename+'_postsolve_display.txt', 'w') as f:
    with contextlib.redirect_stdout(f):
        m.display()

report_infeasibilities(m, "temps/postsolve_violations.log")

pyo.assert_optimal_termination(res)

print(f'capex = {pyo.value(m.fs.capex)}')
print(f'opex = {pyo.value(m.fs.opex)}')
print(f'raw_mats = {pyo.value(m.fs.raw_mats)}')
print(f'revenue = {pyo.value(m.fs.revenue)}')

#visualize flowsheet
m.fs.visualize(filename)

from co2_eor.utility.util_funcs import export_flowsheet_to_excel
export_flowsheet_to_excel(m.fs, 'temps/'+filename+'_solve.xlsx')