import pyomo.environ as pyo
import idaes.core
from pyomo.network import Arc
from pyomo.environ import units
import contextlib

import idaes.models.properties.general_helmholtz as idaesHelmholtz
import idaes.models.unit_models.pressure_changer as idaesPressureChanger

from idaes.core import MomentumBalanceType
from co2_eor.splitter import EnergySplittingType
from co2_eor.mixer import MomentumMixingType

from co2_eor.SSLW import SSLWCosting, SSLWCostingData
from co2_eor.SSLW import CompressorType, CompressorDriveType, CompressorMaterial
from co2_eor.SSLW import VesselMaterial

from co2_eor import pipeline, mixer,splitter, wellpad
from co2_eor.util_funcs import export_compressor_df

from idaes.models.properties.modular_properties.base.generic_property import GenericParameterBlock
from co2_eor.liq_pipe import liqPipe
from co2_eor.wellpattern_r3 import wellpattern
from co2_eor.thermo_config import configuration
from co2_eor.gas_pipe import gasPipe
from co2_eor.processing_facility import processingFacility
from co2_eor.mpf_to_helmholtz import mpf_helmholtz_converter


#model, flowsheet, and parameter block
m = pyo.ConcreteModel()
m.fs = idaes.core.FlowsheetBlock(dynamic=False)
m.fs.props1 = idaesHelmholtz.HelmholtzParameterBlock(
        pure_component="CO2",amount_basis=idaesHelmholtz.AmountBasis.MASS,
        state_vars=idaesHelmholtz.StateVars.PH,phase_presentation=idaesHelmholtz.PhaseType.L,
        #has_phase_equilibrium=False,
        )
m.fs.props2 = GenericParameterBlock(**configuration)

m.fs.costing = SSLWCosting()

compressor_config = {
        "property_package": m.fs.props1,
        "dynamic":False,
        "compressor":True,
        "thermodynamic_assumption":idaesPressureChanger.ThermodynamicAssumption.isentropic,
        }

liqPipe_config = {
        "property_package":m.fs.props1,
        "alpha":5,
        "ambient_temperature":293.15,
        "average_pressure_type":'linear',
        "heat_balance_type":'nonisothermal',
        "average_pressure_weight":0.5,
        "average_temperature_weight":0.0,
        "height_change":0,
        }

mixer_config = {
        "momentum_mixing_type":MomentumMixingType.inequality
        }

splitter_config = {
        "property_package":m.fs.props1,
        "ideal_separation":False,
        "energy_split_basis":EnergySplittingType.equal_molar_enthalpy,
        "momentum_balance_type":MomentumBalanceType.none,
        }

wellpattern_config = {
    "primary_property_package":m.fs.props1,
    "secondary_property_package":m.fs.props2,
    "temperature":300,
    "pressure":413*100000, #413
    "Boi":1.3,
    "GOR":4.42,
    "PC_A":0.4,
    "PC_B":0.6,
    "GB_A":0.62,
    "GB_B":1.12,
    "kovr":2.93E-12,
    "use_correction_factor":False,
    "depth":1000,
    "BHP_max":650*100000, #650?
    "outlet_pressure":30*100000,
    "outlet_temperature":300,
}

gasPipe_config = {
    "property_package":m.fs.props2,
    "average_pressure_type":'nonlinear'
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

#distribution network:
#source branch
m.fs.source_comp = idaesPressureChanger.PressureChanger(**compressor_config)
m.fs.source_comp.costing = idaes.core.UnitModelCostingBlock(**compressor_costing_config)
m.fs.source_comp.efficiency_isentropic.fix(0.85)

m.fs.pipe0 = liqPipe(**liqPipe_config,length=15000)
m.fs.pipe0.costing = idaes.core.UnitModelCostingBlock(**pipe_costing_config)
m.fs.pipe0.diameter.fix(0.8)
m.fs.pipe0.roughness.fix(0.0475e-3)
m.fs.pipe0.max_pressure.fix(650*100000)

m.fs.recycle_comp = idaesPressureChanger.PressureChanger(**compressor_config)
m.fs.recycle_comp.costing = idaes.core.UnitModelCostingBlock(**compressor_costing_config)
m.fs.recycle_comp.efficiency_isentropic.fix(0.85)

m.fs.s0 = Arc(source=m.fs.source_comp.outlet,destination=m.fs.pipe0.inlet)

m.fs.mix0 = mixer(**mixer_config,property_package=m.fs.props1,inlet_list=['from_pipe0','from_recycle_comp'])

m.fs.s01 = Arc(source=m.fs.recycle_comp.outlet,destination=m.fs.mix0.from_recycle_comp)
m.fs.s02 = Arc(source=m.fs.pipe0.outlet,destination=m.fs.mix0.from_pipe0)

m.fs.split0 = splitter(**splitter_config,outlet_list=['to_pipe1','to_pipe3','to_pipe6'])
m.fs.s1 = Arc(source=m.fs.mix0.outlet,destination=m.fs.split0.inlet)

#first branch, wells 1 and 2
m.fs.pipe1 = liqPipe(**liqPipe_config,length=800)
m.fs.pipe1.costing=idaes.core.UnitModelCostingBlock(**pipe_costing_config)
m.fs.pipe1.diameter.fix(0.5)
m.fs.pipe1.roughness.fix(0.0475e-3)
m.fs.pipe1.max_pressure.fix(600*100000)

m.fs.split1 = splitter(**splitter_config,outlet_list=['to_comp1','to_pipe2'])

m.fs.comp1 = idaesPressureChanger.PressureChanger(**compressor_config)
m.fs.comp1.costing = idaes.core.UnitModelCostingBlock(**compressor_costing_config)
m.fs.comp1.efficiency_isentropic.fix(0.85)

m.fs.well1 = wellpattern(**wellpattern_config)
m.fs.well1.HCPV.fix(1.5)

m.fs.pipe2 = liqPipe(**liqPipe_config,length=1600)
m.fs.pipe2.costing=idaes.core.UnitModelCostingBlock(**pipe_costing_config)
m.fs.pipe2.diameter.fix(0.5)
m.fs.pipe2.roughness.fix(0.0475e-3)
m.fs.pipe2.max_pressure.fix(600*100000)

m.fs.comp2 = idaesPressureChanger.PressureChanger(**compressor_config)
m.fs.comp2.costing = idaes.core.UnitModelCostingBlock(**compressor_costing_config)
m.fs.comp2.efficiency_isentropic.fix(0.85)

m.fs.well2 = wellpattern(**wellpattern_config)
m.fs.well2.HCPV.fix(1.25)

m.fs.s2 = Arc(source=m.fs.split0.to_pipe1,destination=m.fs.pipe1.inlet)
m.fs.s3 = Arc(source=m.fs.pipe1.outlet,destination=m.fs.split1.inlet)
m.fs.s4 = Arc(source=m.fs.split1.to_comp1,destination=m.fs.comp1.inlet)
m.fs.s5 = Arc(source=m.fs.comp1.outlet,destination=m.fs.well1.inlet)
m.fs.s6 = Arc(source=m.fs.split1.to_pipe2,destination=m.fs.pipe2.inlet)
m.fs.s7 = Arc(source=m.fs.pipe2.outlet,destination=m.fs.comp2.inlet)
m.fs.s8 = Arc(source=m.fs.comp2.outlet,destination=m.fs.well2.inlet)

#second branch, wells 3 and 4
m.fs.pipe3 = liqPipe(**liqPipe_config,length=3200)
m.fs.pipe3.costing=idaes.core.UnitModelCostingBlock(**pipe_costing_config)
m.fs.pipe3.diameter.fix(0.5)
m.fs.pipe3.roughness.fix(0.0475e-3)
m.fs.pipe3.max_pressure.fix(600*100000)

m.fs.split2 = splitter(**splitter_config,outlet_list=['to_pipe4','to_pipe5'])

m.fs.pipe4 = liqPipe(**liqPipe_config,length=4000)
m.fs.pipe4.costing=idaes.core.UnitModelCostingBlock(**pipe_costing_config)
m.fs.pipe4.diameter.fix(0.5)
m.fs.pipe4.roughness.fix(0.0475e-3)
m.fs.pipe4.max_pressure.fix(600*100000)

m.fs.comp3 = idaesPressureChanger.PressureChanger(**compressor_config)
m.fs.comp3.costing = idaes.core.UnitModelCostingBlock(**compressor_costing_config)
m.fs.comp3.efficiency_isentropic.fix(0.85)

m.fs.well3 = wellpattern(**wellpattern_config)
m.fs.well3.HCPV.fix(1)

m.fs.pipe5 = liqPipe(**liqPipe_config,length=2400)
m.fs.pipe5.costing=idaes.core.UnitModelCostingBlock(**pipe_costing_config)
m.fs.pipe5.diameter.fix(0.5)
m.fs.pipe5.roughness.fix(0.0475e-3)
m.fs.pipe5.max_pressure.fix(600*100000)

m.fs.comp4 = idaesPressureChanger.PressureChanger(**compressor_config)
m.fs.comp4.costing = idaes.core.UnitModelCostingBlock(**compressor_costing_config)
m.fs.comp4.efficiency_isentropic.fix(0.85)

m.fs.well4 = wellpattern(**wellpattern_config)
m.fs.well4.HCPV.fix(0.75)

m.fs.s9 = Arc(source=m.fs.split0.to_pipe3,destination=m.fs.pipe3.inlet)
m.fs.s10 = Arc(source=m.fs.pipe3.outlet,destination=m.fs.split2.inlet)
m.fs.s11 = Arc(source=m.fs.split2.to_pipe4,destination=m.fs.pipe4.inlet)
m.fs.s12 = Arc(source=m.fs.pipe4.outlet,destination=m.fs.comp3.inlet)
m.fs.s13 = Arc(source=m.fs.comp3.outlet,destination=m.fs.well3.inlet)
m.fs.s14 = Arc(source=m.fs.split2.to_pipe5,destination=m.fs.pipe5.inlet)
m.fs.s15 = Arc(source=m.fs.pipe5.outlet,destination=m.fs.comp4.inlet)
m.fs.s16 = Arc(source=m.fs.comp4.outlet,destination=m.fs.well4.inlet)

#third branch, wells 5 and 6
m.fs.pipe6 = liqPipe(**liqPipe_config,length=3200)
m.fs.pipe6.costing=idaes.core.UnitModelCostingBlock(**pipe_costing_config)
m.fs.pipe6.diameter.fix(0.5)
m.fs.pipe6.roughness.fix(0.0475e-3)
m.fs.pipe6.max_pressure.fix(600*100000)

m.fs.split3 = splitter(**splitter_config,outlet_list=['to_pipe7','to_pipe8'])

m.fs.pipe7 = liqPipe(**liqPipe_config,length=2400)
m.fs.pipe7.costing=idaes.core.UnitModelCostingBlock(**pipe_costing_config)
m.fs.pipe7.diameter.fix(0.5)
m.fs.pipe7.roughness.fix(0.0475e-3)
m.fs.pipe7.max_pressure.fix(600*100000)

m.fs.comp5 = idaesPressureChanger.PressureChanger(**compressor_config)
m.fs.comp5.costing = idaes.core.UnitModelCostingBlock(**compressor_costing_config)
m.fs.comp5.efficiency_isentropic.fix(0.85)

m.fs.well5 = wellpattern(**wellpattern_config)
m.fs.well5.HCPV.fix(0.5)

m.fs.pipe8 = liqPipe(**liqPipe_config,length=4000)
m.fs.pipe8.costing=idaes.core.UnitModelCostingBlock(**pipe_costing_config)
m.fs.pipe8.diameter.fix(0.5)
m.fs.pipe8.roughness.fix(0.0475e-3)
m.fs.pipe8.max_pressure.fix(600*100000)

m.fs.comp6 = idaesPressureChanger.PressureChanger(**compressor_config)
m.fs.comp6.costing = idaes.core.UnitModelCostingBlock(**compressor_costing_config)
m.fs.comp6.efficiency_isentropic.fix(0.85)

m.fs.well6 = wellpattern(**wellpattern_config)
m.fs.well6.HCPV.fix(0.25)

m.fs.s17 = Arc(source=m.fs.split0.to_pipe6,destination=m.fs.pipe6.inlet)
m.fs.s18 = Arc(source=m.fs.pipe6.outlet,destination=m.fs.split3.inlet)
m.fs.s19 = Arc(source=m.fs.split3.to_pipe7,destination=m.fs.pipe7.inlet)
m.fs.s20 = Arc(source=m.fs.pipe7.outlet,destination=m.fs.comp5.inlet)
m.fs.s21 = Arc(source=m.fs.comp5.outlet,destination=m.fs.well5.inlet)
m.fs.s22 = Arc(source=m.fs.split3.to_pipe8,destination=m.fs.pipe8.inlet)
m.fs.s23 = Arc(source=m.fs.pipe8.outlet,destination=m.fs.comp6.inlet)
m.fs.s24 = Arc(source=m.fs.comp6.outlet,destination=m.fs.well6.inlet)

#collection network
m.fs.mixPF = mixer(**mixer_config,property_package=m.fs.props2,inlet_list=['from_pipe9','from_pipe12','from_pipe14'])

#collection from wells 1,2, 4
m.fs.pipe9 = gasPipe(**gasPipe_config,length=1000)
m.fs.pipe9.costing=idaes.core.UnitModelCostingBlock(**pipe_costing_config)
m.fs.pipe9.diameter.fix(0.5)
m.fs.pipe9.roughness.fix(0.0475e-3)
m.fs.pipe9.max_pressure.fix(600*100000)

m.fs.mix1 = mixer(**mixer_config,property_package=m.fs.props2,inlet_list=['from_well1','from_pipe10'])

m.fs.pipe10 = gasPipe(**gasPipe_config,length=1600)
m.fs.pipe10.costing=idaes.core.UnitModelCostingBlock(**pipe_costing_config)
m.fs.pipe10.diameter.fix(0.5)
m.fs.pipe10.roughness.fix(0.0475e-3)
m.fs.pipe10.max_pressure.fix(600*100000)

m.fs.mix2 = mixer(**mixer_config,property_package=m.fs.props2,inlet_list=['from_well2','from_pipe11'])

m.fs.pipe11 = gasPipe(**gasPipe_config,length=3200)
m.fs.pipe11.costing=idaes.core.UnitModelCostingBlock(**pipe_costing_config)
m.fs.pipe11.diameter.fix(0.5)
m.fs.pipe11.roughness.fix(0.0475e-3)
m.fs.pipe11.max_pressure.fix(600*100000)

m.fs.s25 = Arc(source=m.fs.well4.outlet,destination=m.fs.pipe11.inlet)
m.fs.s26 = Arc(source=m.fs.pipe11.outlet,destination=m.fs.mix2.from_pipe11)
m.fs.s27 = Arc(source=m.fs.well2.outlet,destination=m.fs.mix2.from_well2)
m.fs.s28 = Arc(source=m.fs.mix2.outlet,destination=m.fs.pipe10.inlet)
m.fs.s29 = Arc(source=m.fs.pipe10.outlet,destination=m.fs.mix1.from_pipe10)
m.fs.s30 = Arc(source=m.fs.well1.outlet,destination=m.fs.mix1.from_well1)
m.fs.s31 = Arc(source=m.fs.mix1.outlet,destination=m.fs.pipe9.inlet)
m.fs.s32 = Arc(source=m.fs.pipe9.outlet,destination=m.fs.mixPF.from_pipe9)

#collection from wells 5,6
m.fs.pipe12 = gasPipe(**gasPipe_config,length=2000)
m.fs.pipe12.costing=idaes.core.UnitModelCostingBlock(**pipe_costing_config)
m.fs.pipe12.diameter.fix(0.5)
m.fs.pipe12.roughness.fix(0.0475e-3)
m.fs.pipe12.max_pressure.fix(600*100000)

m.fs.mix3 = mixer(**mixer_config,property_package=m.fs.props2,inlet_list=['from_well6','from_pipe13'])

m.fs.pipe13 = gasPipe(**gasPipe_config,length=5500)
m.fs.pipe13.costing=idaes.core.UnitModelCostingBlock(**pipe_costing_config)
m.fs.pipe13.diameter.fix(0.5)
m.fs.pipe13.roughness.fix(0.0475e-3)
m.fs.pipe13.max_pressure.fix(600*100000)

m.fs.s33 = Arc(source=m.fs.well5.outlet,destination=m.fs.pipe13.inlet)
m.fs.s34 = Arc(source=m.fs.pipe13.outlet,destination=m.fs.mix3.from_pipe13)
m.fs.s35 = Arc(source=m.fs.well6.outlet,destination=m.fs.mix3.from_well6)
m.fs.s36 = Arc(source=m.fs.mix3.outlet,destination=m.fs.pipe12.inlet)
m.fs.s37 = Arc(source=m.fs.pipe12.outlet,destination=m.fs.mixPF.from_pipe12)

#collection from well 3
m.fs.pipe14 = gasPipe(**gasPipe_config,length=2800)
m.fs.pipe14.costing=idaes.core.UnitModelCostingBlock(**pipe_costing_config)
m.fs.pipe14.diameter.fix(0.5)
m.fs.pipe14.roughness.fix(0.0475e-3)
m.fs.pipe14.max_pressure.fix(600*100000)

m.fs.s38 = Arc(source=m.fs.well3.outlet,destination=m.fs.pipe14.inlet)
m.fs.s39 = Arc(source=m.fs.pipe14.outlet,destination=m.fs.mixPF.from_pipe14)

#processing facility to recycle compressor
m.fs.procFac = processingFacility(
    property_package=m.fs.props2,
    has_phase_equilibrium=False,
    key_component='co2',
    minimum_inlet_pressure=5*100000,
    recycle_pressure=2*100000,
)

m.fs.mpfHelmholtzConverter = mpf_helmholtz_converter(
    property_package_in = m.fs.props2,
    property_package_out = m.fs.props1,
    conversion_type = 'convert_co2'
)

m.fs.s40 = Arc(source=m.fs.mixPF.outlet,destination=m.fs.procFac.inlet)
m.fs.s41 = Arc(source=m.fs.procFac.recycle,destination=m.fs.mpfHelmholtzConverter.inlet)
m.fs.s42 = Arc(source=m.fs.mpfHelmholtzConverter.outlet,destination=m.fs.recycle_comp.inlet)

#expand arcs
pyo.TransformationFactory("network.expand_arcs").apply_to(m)

#check degrees of freedom
print(f'D.o.F before specifying anything={idaes.core.util.model_statistics.degrees_of_freedom(m)}')

#inlet specifications
m.fs.source_comp.inlet.pressure[0].fix(100*100000)
m.fs.source_comp.inlet.enth_mass[0].fix(m.fs.props1.htpx(T=298*units.K,p=100*100000*units.Pa,amount_basis=idaesHelmholtz.AmountBasis.MASS))

#initialization dof

m.fs.source_comp.outlet.pressure[0].fix(600*100000)
m.fs.source_comp.inlet.flow_mass[0].fix(1)

m.fs.recycle_comp.outlet.pressure[0].fix(600*100000)

m.fs.comp1.outlet.pressure[0].fix(600*100000)
m.fs.comp2.outlet.pressure[0].fix(600*100000)
m.fs.comp3.outlet.pressure[0].fix(600*100000)
m.fs.comp4.outlet.pressure[0].fix(600*100000)
m.fs.comp5.outlet.pressure[0].fix(600*100000)
m.fs.comp6.outlet.pressure[0].fix(600*100000)

with open('temps/bakken_network_pprint.txt', 'w') as f:
    with contextlib.redirect_stdout(f):
        m.pprint()

#initialize

from pyomo.network import SequentialDecomposition
from idaes.core.util.initialization import propagate_state
seq = SequentialDecomposition()
seq.options.select_tear_method = "heuristic"
seq.options.tear_method = "Wegstein"
seq.options.iterLim = 5

G = seq.create_graph(m)
heauristic_tear_set = seq.tear_set_arcs(G,method="heuristic")
order = seq.calculation_order(G)

print("tear set")
for o in heauristic_tear_set:
    print(o.name)
print("calculation order")
for o in order:
    print(o[0].name)

tear_guesses = {
    "flow_mass":{
        (0):1,
    },
    "temperature":{0:300},
    "pressure":{0:600*100000},
}
seq.set_guesses_for(m.fs.mix0.from_recycle_comp,tear_guesses)

from co2_eor.util_funcs import fs_initializer_function

seq.run(m,fs_initializer_function)
print(f'initialization done')

#visualize flowsheet
m.fs.visualize("bakken_network")
input()

