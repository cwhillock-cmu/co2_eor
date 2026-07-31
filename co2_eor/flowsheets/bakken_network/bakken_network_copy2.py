import pyomo.environ as pyo
import idaes.core
from pyomo.network import Arc
from pyomo.environ import units
import contextlib
import numpy as np

import idaes.models.properties.general_helmholtz as idaesHelmholtz
import idaes.models.unit_models.pressure_changer as idaesPressureChanger

from idaes.core import MomentumBalanceType
from co2_eor.splitter_unit import EnergySplittingType
from co2_eor.mixer_unit import MomentumMixingType

from co2_eor.SSLW import SSLWCosting, SSLWCostingData
from co2_eor.SSLW import CompressorType, CompressorDriveType, CompressorMaterial
from co2_eor.SSLW import VesselMaterial

from co2_eor import pipeline, mixer,splitter, wellpad, heater
from co2_eor.util_funcs import export_compressor_df

from idaes.models.properties.modular_properties.base.generic_property import GenericParameterBlock
from co2_eor.liq_pipe import liqPipe
from co2_eor.wellpattern_HH import wellpattern
from co2_eor.MPF.thermo_config import configuration_vap
from co2_eor.gas_pipe import gasPipe
from co2_eor.processing_facility import processingFacility
from co2_eor.mpf_to_helmholtz import mpf_helmholtz_converter

from co2_eor.util_funcs import ipopt, conopt


#model, flowsheet, and parameter block
m = pyo.ConcreteModel()
m.fs = idaes.core.FlowsheetBlock(dynamic=False)
m.fs.props1 = idaesHelmholtz.HelmholtzParameterBlock(
        pure_component="CO2",amount_basis=idaesHelmholtz.AmountBasis.MOLE,
        state_vars=idaesHelmholtz.StateVars.PH,phase_presentation=idaesHelmholtz.PhaseType.MIX,
        #has_phase_equilibrium=False,
        )
m.fs.props2 = GenericParameterBlock(**configuration_vap)

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
        "allowable_stress":2000*100000,
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
    "kovr":2.93E-14,
    "use_correction_factor":False,
    "depth":1000,
    "BHP_max":650*100000, #650?
    "outlet_pressure":30*100000,
    "outlet_temperature":300,
}

heater_config = {
    "property_package":m.fs.props1,
}

gasPipe_config = {
    "property_package":m.fs.props2,
    "average_pressure_type":'nonlinear',
    "allowable_stress":2000*100000,
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

m.fs.chiller = heater(**heater_config)
m.fs.chiller.outlet_temp_eq = pyo.Constraint(expr=m.fs.chiller.control_volume.properties_out[0].temperature<=320)
m.fs.chiller.outlet_temp_eq2 = pyo.Constraint(expr=m.fs.chiller.control_volume.properties_out[0].temperature>=290)

m.fs.s_source_comp_pipe0 = Arc(source=m.fs.source_comp.outlet,destination=m.fs.pipe0.inlet)

m.fs.mix0 = mixer(**mixer_config,property_package=m.fs.props1,inlet_list=['from_pipe0','from_chiller'])

m.fs.s_recycle_comp_chiller = Arc(source=m.fs.recycle_comp.outlet,destination=m.fs.chiller.inlet)
m.fs.s_chiller_mix0 = Arc(source=m.fs.chiller.outlet,destination=m.fs.mix0.from_chiller)
m.fs.s_pipe0_mix0 = Arc(source=m.fs.pipe0.outlet,destination=m.fs.mix0.from_pipe0)

m.fs.split0 = splitter(**splitter_config,outlet_list=['to_pipe1','to_pipe3','to_pipe6'])
m.fs.s_mix0_split0 = Arc(source=m.fs.mix0.outlet,destination=m.fs.split0.inlet)

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
m.fs.comp1.pressure_minimum_eq = pyo.Constraint(expr=m.fs.comp1.inlet.pressure[0]>=80*100000)

m.fs.well1 = wellpattern(**wellpattern_config)
m.fs.well1.HCPV.fix(1.1)

m.fs.pipe2 = liqPipe(**liqPipe_config,length=1600)
m.fs.pipe2.costing=idaes.core.UnitModelCostingBlock(**pipe_costing_config)
m.fs.pipe2.diameter.fix(0.5)
m.fs.pipe2.roughness.fix(0.0475e-3)
m.fs.pipe2.max_pressure.fix(600*100000)

m.fs.comp2 = idaesPressureChanger.PressureChanger(**compressor_config)
m.fs.comp2.costing = idaes.core.UnitModelCostingBlock(**compressor_costing_config)
m.fs.comp2.efficiency_isentropic.fix(0.85)
m.fs.comp2.pressure_minimum_eq = pyo.Constraint(expr=m.fs.comp2.inlet.pressure[0]>=80*100000)

m.fs.well2 = wellpattern(**wellpattern_config)
m.fs.well2.HCPV.fix(1.0)

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
m.fs.comp3.pressure_minimum_eq = pyo.Constraint(expr=m.fs.comp3.inlet.pressure[0]>=80*100000)

m.fs.well3 = wellpattern(**wellpattern_config)
m.fs.well3.HCPV.fix(0.8)

m.fs.pipe5 = liqPipe(**liqPipe_config,length=2400)
m.fs.pipe5.costing=idaes.core.UnitModelCostingBlock(**pipe_costing_config)
m.fs.pipe5.diameter.fix(0.5)
m.fs.pipe5.roughness.fix(0.0475e-3)
m.fs.pipe5.max_pressure.fix(600*100000)

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
m.fs.comp5.pressure_minimum_eq = pyo.Constraint(expr=m.fs.comp5.inlet.pressure[0]>=80*100000)

m.fs.well5 = wellpattern(**wellpattern_config)
m.fs.well5.HCPV.fix(0.9)

m.fs.pipe8 = liqPipe(**liqPipe_config,length=4000)
m.fs.pipe8.costing=idaes.core.UnitModelCostingBlock(**pipe_costing_config)
m.fs.pipe8.diameter.fix(1)
m.fs.pipe8.roughness.fix(0.0475e-3)
m.fs.pipe8.max_pressure.fix(600*100000)

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
m.fs.pipe12.diameter.fix(0.5)
m.fs.pipe12.roughness.fix(0.0475e-3)
m.fs.pipe12.max_pressure.fix(600*100000)

m.fs.mix3 = mixer(**mixer_config,property_package=m.fs.props2,inlet_list=['from_well6','from_pipe13'])

m.fs.pipe13 = gasPipe(**gasPipe_config,length=5500)
m.fs.pipe13.costing=idaes.core.UnitModelCostingBlock(**pipe_costing_config)
m.fs.pipe13.diameter.fix(0.5)
m.fs.pipe13.roughness.fix(0.0475e-3)
m.fs.pipe13.max_pressure.fix(600*100000)

m.fs.s_well5_pipe13 = Arc(source=m.fs.well5.outlet,destination=m.fs.pipe13.inlet)
m.fs.s_pipe13_mix3 = Arc(source=m.fs.pipe13.outlet,destination=m.fs.mix3.from_pipe13)
m.fs.s_well6_mix3 = Arc(source=m.fs.well6.outlet,destination=m.fs.mix3.from_well6)
m.fs.s_mix3_pipe12 = Arc(source=m.fs.mix3.outlet,destination=m.fs.pipe12.inlet)
m.fs.s_pipe12_mixPF = Arc(source=m.fs.pipe12.outlet,destination=m.fs.mixPF.from_pipe12)

#collection from well 3
m.fs.pipe14 = gasPipe(**gasPipe_config,length=2800)
m.fs.pipe14.costing=idaes.core.UnitModelCostingBlock(**pipe_costing_config)
m.fs.pipe14.diameter.fix(0.5)
m.fs.pipe14.roughness.fix(0.0475e-3)
m.fs.pipe14.max_pressure.fix(600*100000)

m.fs.s_well3_pipe14 = Arc(source=m.fs.well3.outlet,destination=m.fs.pipe14.inlet)
m.fs.s_pipe14_mixPF = Arc(source=m.fs.pipe14.outlet,destination=m.fs.mixPF.from_pipe14)

#processing facility to recycle compressor
m.fs.procFac = processingFacility(
    property_package=m.fs.props2,
    has_phase_equilibrium=False,
    key_component='co2',
    minimum_inlet_pressure=5*100000,
    recycle_pressure=5*100000,
)
"""
m.fs.procFac.costing = customCostingBlock(costing_method=cost_processing_facility_simple, costing_method_args={
    "F_in":m.fs.procFac.control_volume.properties_in[0].flow_mass_comp['co2'],
    "F_recycle":m.fs.procFac.control_volume.properties_recycle[0].flow_mass_comp['co2'],
    "oil_price":70,
})
"""

m.fs.mpfHelmholtzConverter = mpf_helmholtz_converter(
    property_package_in = m.fs.props2,
    property_package_out = m.fs.props1,
    conversion_type = 'convert_all_mol'
)

m.fs.purge_splitter = splitter(**splitter_config,outlet_list=['purge','to_recycle_comp'])
m.fs.purge_splitter.positive_purge_eq = pyo.Constraint(expr=m.fs.purge_splitter.purge.flow_mol[0]>=0) #enforce material does not go into system from purge stream

m.fs.s_mixPF_procFac = Arc(source=m.fs.mixPF.outlet,destination=m.fs.procFac.inlet)
m.fs.s_procFac_eos_converter = Arc(source=m.fs.procFac.recycle,destination=m.fs.mpfHelmholtzConverter.inlet)
m.fs.s_eos_converter_purge_splitter = Arc(source=m.fs.mpfHelmholtzConverter.outlet,destination=m.fs.purge_splitter.inlet)
m.fs.s_purge_splitter_recycle_comp = Arc(source=m.fs.purge_splitter.to_recycle_comp,destination=m.fs.recycle_comp.inlet)

#expand arcs
pyo.TransformationFactory("network.expand_arcs").apply_to(m)

#check degrees of freedom
print(f'D.o.F before specifying anything={idaes.core.util.model_statistics.degrees_of_freedom(m)}')

#inlet specifications
m.fs.source_comp.inlet.pressure[0].fix(100*100000)
m.fs.source_comp.inlet.enth_mol[0].fix(m.fs.props1.htpx(T=298*units.K,p=100*100000*units.Pa,amount_basis=idaesHelmholtz.AmountBasis.MOLE))

#initialization dof

m.fs.source_comp.outlet.pressure[0].fix(200*100000)
m.fs.source_comp.inlet.flow_mol[0].fix(80)

m.fs.recycle_comp.outlet.pressure[0].fix(200*100000)

for i in range(1,7):
    getattr(m.fs, "comp"+str(i)).outlet.pressure[0].fix(650*100000)

m.fs.purge_splitter.split_fraction[0,'purge'].fix(1e-14)

with open('temps/bakken_network_preinitialization1_pprint.txt', 'w') as f:
    with contextlib.redirect_stdout(f):
        m.pprint()

#initialization 1

from pyomo.network import SequentialDecomposition
from idaes.core.util.initialization import propagate_state
seq = SequentialDecomposition()
seq.options.select_tear_method = "heuristic"
seq.options.tear_method = "Wegstein"
seq.options.iterLim = 10

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
    "flow_mol":{
        (0):600,
    },
    "enth_mol":{0:m.fs.props1.htpx(T=300*units.K,p=200*100000*units.Pa,amount_basis=idaesHelmholtz.AmountBasis.MOLE)},
    "pressure":{0:200*100000},
}
seq.set_guesses_for(m.fs.chiller.inlet,tear_guesses)

from co2_eor.util_funcs import custom_initialize_unit
def fs_custom_SD_initializer(unit):
    print(f'initializing {unit.name}')
    comp_list = []
    well_list = []
    for i in range(1,7):
        comp_list.append("fs.comp"+str(i))
        well_list.append("fs.well"+str(i))
    if unit.name in comp_list:
        number = unit.name[-1]
        m.temp_fs = pyo.Block()
        m.temp_fs.comp = pyo.Reference(getattr(m.fs,"comp"+number))
        m.temp_fs.well = pyo.Reference(getattr(m.fs,"well"+number))
        m.temp_fs.stream = pyo.Reference(getattr(m.fs,"s"+"_comp"+number+"_well"+number))
        m.temp_fs.stream_expanded = pyo.Reference(getattr(m.fs,"s"+"_comp"+number+"_well"+number+"_expanded"))
        m.temp_fs.well[None].activate_feasibility_problem()
        ipopt.options['nlp_scaling_method'] = 'user-scaling'
        try:
            res=ipopt.solve(m.temp_fs,tee=False)
        except ValueError:
            with open('temps/temp_fs_pprint.txt', 'w') as f:
                with contextlib.redirect_stdout(f):
                    m.temp_fs.pprint()
            input('paused')
            pass
        ipopt.options['nlp_scaling_method'] = 'gradient-based'
        m.temp_fs.well[None].deactivate_feasibility_problem()
        m.del_component(m.temp_fs)

    elif unit.name in well_list:
        return
    else:
        custom_initialize_unit(unit)

seq.run(m,fs_custom_SD_initializer)
print(f'initialization 1 done')

#unfix initialization 1 degrees of freedom
m.fs.source_comp.inlet.flow_mol[0].unfix()
m.fs.source_comp.outlet.pressure[0].unfix()
m.fs.recycle_comp.outlet.pressure[0].unfix()
#m.fs.purge_splitter.split_fraction[0,'purge'].unfix()
#for i in range(1,7):
#    getattr(m.fs, "comp"+str(i)).outlet.pressure[0].unfix()
#for i in range(1,15):
#    getattr(m.fs, "pipe"+str(i)).diameter.unfix()

#create production target
m.fs.prod_target = pyo.Param(mutable=True,initialize=570000) #STB/yr
#production target constraint
@m.fs.Constraint()
def production_target_constraint(fs):
    rhs = 0
    for i in range(1,7):
        rhs += getattr(fs,"well"+str(i)).q_OIL_PROD
    return fs.prod_target / 365 <= rhs * 3600 * 24
#m.fs.production_target_constraint.deactivate()

#initialization 2

for i in range(1,15):
    getattr(m.fs, "pipe"+str(i)).activate_slack_variables()
#for i in range(1,7):
#    getattr(m.fs,"well"+str(i)).activate_slack_variables()

@m.fs.Objective()
def flowsheet_feasibility_objective(fs):
    expr = 0
    for i in range(1,15):
        expr += getattr(fs, "pipe"+str(i)).feasibility_expression
    #for i in range(1,7):
    #    expr += getattr(fs,"well"+str(i)).feasibility_expression
    return expr

#unfix initialization 2 dof
m.fs.source_comp.inlet.flow_mol[0].unfix()
m.fs.source_comp.outlet.pressure[0].unfix()
m.fs.recycle_comp.outlet.pressure[0].unfix()

with open('temps/bakken_network_preinitialization2_pprint.txt', 'w') as f:
    with contextlib.redirect_stdout(f):
        m.pprint()

#scale model
scaled_m = pyo.TransformationFactory("core.scale_model").create_using(m)
#solve flowsheet
ipopt.options['acceptable_tol']=1E-4
ipopt.options['tol']=1E-6
ipopt.options['max_iter'] = 3000
ipopt.options['linear_solver']='ma97'
try: 
    res=ipopt.solve(scaled_m,tee=False)
    pass
except ValueError:
    pass
ipopt.options['acceptable_tol']=1E-6
ipopt.options['tol']=1E-8
ipopt.options['linear_solver']='ma27'
#unscale model
pyo.TransformationFactory("core.scale_model").propagate_solution(scaled_m,m)

#deactivate flowsheet feasibility problem
m.fs.flowsheet_feasibility_objective.deactivate()
for i in range(1,15):
    getattr(m.fs, "pipe"+str(i)).deactivate_feasibility_problem()
for i in range(1,7):
    getattr(m.fs,"well"+str(i)).deactivate_feasibility_problem()

print(f'initialization 2 done')

with open('temps/bakken_network_postinitialization2_display.txt', 'w') as f:
    with contextlib.redirect_stdout(f):
        m.display()

#initialization 3
m.fs.source_comp.inlet.flow_mol[0].fix()
m.fs.source_comp.outlet.pressure[0].fix()
m.fs.recycle_comp.outlet.pressure[0].fix()
m.fs.purge_splitter.split_fraction[0,'purge'].fix()

def SD_solve(unit):
    print(f'solving {unit.name}')
    ipopt.options['acceptable_tol']=1E-4
    ipopt.options['tol']=1E-6
    ipopt.options['max_iter'] = 300
    comp_list = []
    well_list = []
    for i in range(1,7):
        comp_list.append("fs.comp"+str(i))
        well_list.append("fs.well"+str(i))
    if unit.name in comp_list:
        number = unit.name[-1]
        m.temp_fs = pyo.Block()
        m.temp_fs.comp = pyo.Reference(getattr(m.fs,"comp"+number))
        m.temp_fs.well = pyo.Reference(getattr(m.fs,"well"+number))
        m.temp_fs.stream = pyo.Reference(getattr(m.fs,"s"+"_comp"+number+"_well"+number))
        m.temp_fs.stream_expanded = pyo.Reference(getattr(m.fs,"s"+"_comp"+number+"_well"+number+"_expanded"))
        ipopt.options['nlp_scaling_method'] = 'user-scaling'
        try:
            res=ipopt.solve(m.temp_fs,tee=False)
        except ValueError:
            with open('temps/temp_fs_pprint.txt', 'w') as f:
                with contextlib.redirect_stdout(f):
                    m.temp_fs.pprint()
            #input('paused')
            pass
        ipopt.options['nlp_scaling_method'] = 'gradient-based'
        m.del_component(m.temp_fs)

    elif unit.name in well_list:
        return
    else:
        scaled_unit = pyo.TransformationFactory("core.scale_model").create_using(unit)
        if unit.name[:6] == 'fs.mix':
            ipopt.options['max_iter'] = 9
        try:
            res=ipopt.solve(scaled_unit,tee=False)
            pyo.assert_optimal_termination(res)
        except (ValueError, RuntimeError) as error:
            if unit.name[:6] != 'fs.mix':
                unit.initialize()
            #unit.initialize()
            #with open('temps/SD_solve_unit_pprint.txt', 'w') as f:
            #    with contextlib.redirect_stdout(f):
            #        unit.pprint()
            #input('paused')
            pass
        pyo.TransformationFactory("core.scale_model").propagate_solution(scaled_unit,unit)
    ipopt.options['max_iter'] = 3000
    ipopt.options['acceptable_tol']=1E-6
    ipopt.options['tol']=1E-8

seq.options.iterLim = 10
seq.run(m,SD_solve)
print(f'initialization 3 done')

with open('temps/bakken_network_postinitialization3_display.txt', 'w') as f:
    with contextlib.redirect_stdout(f):
        m.display()

#unfix all initialization degrees of freedom
m.fs.source_comp.inlet.flow_mol[0].unfix()
m.fs.source_comp.outlet.pressure[0].unfix()
m.fs.recycle_comp.outlet.pressure[0].unfix()
m.fs.purge_splitter.split_fraction[0,'purge'].unfix()
for i in range(1,7):
    getattr(m.fs, "comp"+str(i)).outlet.pressure[0].unfix()
for i in range(0,15):
    getattr(m.fs, "pipe"+str(i)).diameter.unfix()

#create objective function terms
@m.fs.Expression()
def opex(fs):
    expr = fs.source_comp.work_mechanical[0] + fs.recycle_comp.work_mechanical[0] - 0.2*fs.chiller.heat_duty[0]
    for i in range(1,7):
        expr += getattr(fs,"comp"+str(i)).work_mechanical[0]
    return 365*24*3600* 4E-8*expr + fs.procFac.costing.opex/3600/24/365

m.fs.raw_mats = pyo.Expression(
    expr=365*24*3600* 0.050*m.fs.pipe0.control_volume.properties_in[0].flow_mass
)

@m.fs.Expression()
def revenue(fs):
    expr = 0
    for i in range(1,7):
        expr += getattr(fs,"well"+str(i)).q_OIL_PROD
    return 365*24*3600* 70*expr

@m.fs.Expression()
def capex(fs):
    expr = fs.source_comp.costing.capital_cost + fs.recycle_comp.costing.capital_cost + fs.procFac.costing.capital_cost
    for i in range(1,7):
        expr += getattr(fs,"comp"+str(i)).costing.capital_cost
    for i in range(0,15):
        expr += getattr(fs,"pipe"+str(i)).costing.capital_cost
    return expr

m.fs.obj_with_revenue = pyo.Objective(
    expr=(m.fs.capex*0.1+m.fs.opex+m.fs.raw_mats-m.fs.revenue)/10000
)
m.fs.obj_no_revenue = pyo.Objective(
    expr=(m.fs.capex*0.1+m.fs.opex+m.fs.raw_mats)/10000
)
m.fs.obj_all_revenue = pyo.Objective(
    expr=(-m.fs.revenue)/10000
)
m.fs.obj_with_revenue.deactivate()
m.fs.obj_all_revenue.deactivate()
m.fs.obj_no_revenue.activate()
m.fs.production_target_constraint.activate()

with open('temps/bakken_network_presolve_pprint.txt', 'w') as f:
    with contextlib.redirect_stdout(f):
        m.pprint()

#scale model
scaled_m = pyo.TransformationFactory("core.scale_model").create_using(m)
#solve flowsheet
ipopt.options['max_iter'] = 4000
ipopt.options['linear_solver']='ma27'
ipopt.options['acceptable_tol']=1E-6
ipopt.options['tol']=1E-8
ipopt.options['nlp_scaling_method']='gradient-based'
ipopt.options['expect_infeasible_problem']='yes'
ipopt.options['OF_evaluate_orig_obj_at_resto_trial']='no'
ipopt.options['OF_start_with_resto']='yes'
#ipopt.options['OF_soft_resto_pderror_reduction_factor']=0
res=ipopt.solve(scaled_m,tee=True)
#unscale model
pyo.TransformationFactory("core.scale_model").propagate_solution(scaled_m,m)

with open('temps/bakken_network_postsolve_display.txt', 'w') as f:
    with contextlib.redirect_stdout(f):
        m.display()

pyo.assert_optimal_termination(res)

print(f'capex = {pyo.value(m.fs.capex)}')
print(f'opex = {pyo.value(m.fs.opex)}')
print(f'raw_mats = {pyo.value(m.fs.raw_mats)}')
print(f'revenue = {pyo.value(m.fs.revenue)}')

#visualize flowsheet
#m.fs.visualize("bakken_network")

from co2_eor.util_funcs import export_flowsheet_to_excel
export_flowsheet_to_excel(m.fs, 'temps/bakken_network_solve.xlsx')

#pause program
#"""
import readchar
print("Press any key to continue...")
key = readchar.readkey() 
print(f"Resumed after pressing: {key}")
#"""

m.fs.obj_with_revenue.deactivate()
#m.fs.obj_no_revenue.activate()
m.fs.obj_all_revenue.deactivate()
m.fs.production_target_constraint.activate()

from pyomo.common.errors import ApplicationError
#prod_targets = np.linspace(1200000,100000,56)
#prod_targets = np.flip(np.arange(50000,675000,25000))
prod_targets = np.flip(np.arange(425000,605000,5000))
ipopt.options['max_iter'] = 2500
ipopt.options['linear_solver']='ma27'
ipopt.options['acceptable_tol']=1E-6
ipopt.options['tol']=1E-8
ipopt.options['nlp_scaling_method']='gradient-based'
ipopt.options['expect_infeasible_problem']='yes'
ipopt.options['OF_evaluate_orig_obj_at_resto_trial']='no'
ipopt.options['OF_start_with_resto']='yes'
#ipopt.options['nlp_scaling_method']='none'
conopt.options['iterlim']=500
for prod_target in prod_targets:
    run_conopt = False
    m.fs.prod_target=prod_target
    print(f'starting solve prod target = {prod_target}')
    #"""
    try:
        scaled_m = pyo.TransformationFactory("core.scale_model").create_using(m)
        with open(rf'temps/bakken_sweep2/bakken_{int(prod_target)}_ipopt_output.txt', 'w') as f:
            with contextlib.redirect_stdout(f):
                res = ipopt.solve(scaled_m,tee=True)
        if res.solver.termination_condition == pyo.TerminationCondition.optimal:
            pyo.TransformationFactory("core.scale_model").propagate_solution(scaled_m,m)
            with open(rf'temps/bakken_sweep2/bakken_{int(prod_target)}_postsolve_display.txt', 'w') as f:
                with contextlib.redirect_stdout(f):
                    m.display()
            print(f'target {prod_target} solved')
            export_flowsheet_to_excel(m.fs, rf'temps/bakken_sweep2/bakken_{int(prod_target)}_export.xlsx')
            print('ipopt converged')
        else:
            run_conopt = True
            print('ipopt not optimal')
    except (ValueError, ApplicationError) as error:
        print(f'ipopt failed with {error}')
        run_conopt = True
        pass
    #"""
    if run_conopt:
        print('switching to conopt')
        scaled_m = pyo.TransformationFactory("core.scale_model").create_using(m)
        with open(rf'temps/bakken_sweep2/bakken_{int(prod_target)}_conopt_output.txt', 'w') as f:
            with contextlib.redirect_stdout(f):
                res = conopt.solve(scaled_m,tee=True)
        if res.solver.termination_condition == pyo.TerminationCondition.optimal:
            pyo.TransformationFactory("core.scale_model").propagate_solution(scaled_m,m)
            with open(rf'temps/bakken_sweep2/bakken_{int(prod_target)}_postsolve_display.txt', 'w') as f:
                with contextlib.redirect_stdout(f):
                    m.display()
            print(f'target {prod_target} solved')
            export_flowsheet_to_excel(m.fs, rf'temps/bakken_sweep2/bakken_{int(prod_target)}_export.xlsx')
            print(f'conopt converged')
        else:
            print(f'target {prod_target} failed')
        
