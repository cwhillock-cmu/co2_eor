import pyomo.environ as pyo
import idaes.core
from pyomo.network import Arc

import idaes.models.properties.general_helmholtz as idaesHelmholtz
import idaes.models.unit_models.pressure_changer as idaesPressureChanger

from idaes.core import MomentumBalanceType
from co2_eor.splitter import EnergySplittingType
from co2_eor.mixer import MomentumMixingType

from idaes.models.costing.SSLW import SSLWCosting, SSLWCostingData
from idaes.models.costing.SSLW import CompressorType, CompressorDriveType, CompressorMaterial
from idaes.models.costing.SSLW import VesselMaterial

from co2_eor import pipeline, wellpad, mixer,splitter
from co2_eor.util_funcs import export_compressor_df

#model, flowsheet, and parameter block
m = pyo.ConcreteModel()
m.fs = idaes.core.FlowsheetBlock(dynamic=False)
m.fs.props = idaesHelmholtz.HelmholtzParameterBlock(
        pure_component="CO2",amount_basis=idaesHelmholtz.AmountBasis.MASS,
        state_vars=idaesHelmholtz.StateVars.PH,phase_presentation=idaesHelmholtz.PhaseType.LG,
        #has_phase_equilibrium=False,
        )

m.fs.costing = SSLWCosting()

compressor_config = {
        "property_package": m.fs.props,
        "dynamic":False,
        "compressor":True,
        "thermodynamic_assumption":idaesPressureChanger.ThermodynamicAssumption.isentropic,
        }

pipeline_config = {
        "property_package":m.fs.props,
        "alpha":5,
        "ambient_temperature":293.15,
        "average_pressure_type":'linear',
        "heat_balance_type":'nonisothermal',
        "average_pressure_weight":0.5,
        "average_temperature_weight":0.0,
        "height_change":0,
        }

mixer_config = {
        "property_package":m.fs.props,
        "momentum_mixing_type":MomentumMixingType.inequality
        }

splitter_config = {
        "property_package":m.fs.props,
        "ideal_separation":False,
        "energy_split_basis":EnergySplittingType.equal_molar_enthalpy,
        "momentum_balance_type":MomentumBalanceType.none,
        }

welltype1_config = {
        "property_package":m.fs.props,
        "temperature":306,
        "pressure":140*100000,
        "Boi":1.2,
        "GOR":0.6,
        "PC_A":0.4,
        "PC_B":0.6,
        "GB_A":0.62,
        "GB_B":1.12,
        "kovr":1.5E-14,
        "use_correction_factor":False,
        "BHP_max":300*100000,
        "multiplier":3,
        "depth":1000,
        }

welltype2_config = {
        "property_package":m.fs.props,
        "temperature":295,
        "pressure":150*100000,
        "Boi":1.3,
        "GOR":0.7,
        "PC_A":0.3,
        "PC_B":0.5,
        "GB_A":0.53,
        "GB_B":1.18,
        "kovr":8.7E-15,
        "use_correction_factor":False,
        "BHP_max":320*100000,
        "multiplier":2,
        "depth":1000,
        }

welltype3_config={
        "property_package":m.fs.props,
        "temperature":340,
        "pressure":120*100000,
        "Boi":1.4,
        "GOR":0.75,
        "PC_A":0.25,
        "PC_B":0.3,
        "GB_A":0.67,
        "GB_B":1.21,
        "kovr":0.2E-14,
        "use_correction_factor":False,
        "BHP_max":340*100000,
        "multiplier":2,
        "depth":1000,
        }

pipeline_costing_config={
    "flowsheet_costing_block":m.fs.costing,
    "costing_method":SSLWCostingData.cost_horizontal_vessel,
    "costing_method_arguments":{
        "material_type":VesselMaterial.StainlessSteel316,
        "include_platforms_ladders":False,
    }
}

#compressor 1
m.fs.comp1 = idaesPressureChanger.PressureChanger(**compressor_config)
m.fs.comp1.efficiency_isentropic.fix(0.85)
m.fs.comp1.costing = idaes.core.UnitModelCostingBlock(
    flowsheet_costing_block=m.fs.costing,
    costing_method=SSLWCostingData.cost_compressor,
    costing_method_arguments={
        "compressor_type":CompressorType.Centrifugal,
        "drive_type":CompressorDriveType.ElectricMotor,
        "material_type":CompressorMaterial.StainlessSteel
    }
)

m.fs.pipe1 = pipeline(**pipeline_config,length=3000)
m.fs.pipe1.costing = idaes.core.UnitModelCostingBlock(**pipeline_costing_config)
m.fs.pipe1.diameter.fix(0.2)
m.fs.pipe1.roughness.fix(0.0475e-3)

m.fs.s1 = Arc(source=m.fs.comp1.outlet,destination=m.fs.pipe1.inlet)

m.fs.split1 = splitter(**splitter_config,outlet_list=['to_pipe2','to_pipe4'])

m.fs.s2 = Arc(source=m.fs.pipe1.outlet,destination=m.fs.split1.inlet)

m.fs.pipe2 = pipeline(**pipeline_config,length=1600)
m.fs.pipe2.costing=idaes.core.UnitModelCostingBlock(**pipeline_costing_config)
m.fs.pipe2.diameter.fix(0.15)
m.fs.pipe2.roughness.fix(0.0475e-3)

m.fs.s3 = Arc(source=m.fs.split1.to_pipe2,destination=m.fs.pipe2.inlet)

m.fs.split2 = splitter(**splitter_config,outlet_list=['to_well1','to_pipe3'])

m.fs.s4 = Arc(source=m.fs.pipe2.outlet,destination=m.fs.split2.inlet)

m.fs.well1 = wellpad(**welltype3_config)
m.fs.well1.HCPV.fix(0.5)

m.fs.s5 = Arc(source=m.fs.split2.to_well1,destination=m.fs.well1.inlet)

m.fs.pipe3 = pipeline(**pipeline_config,length=5600)
m.fs.pipe3.costing = idaes.core.UnitModelCostingBlock(**pipeline_costing_config)
m.fs.pipe3.diameter.fix(0.1)
m.fs.pipe3.roughness.fix(0.0475e-3)

m.fs.s6 = Arc(source=m.fs.split2.to_pipe3,destination=m.fs.pipe3.inlet)

m.fs.well2 = wellpad(**welltype1_config)
m.fs.well2.HCPV.fix(0.5)

m.fs.s7 = Arc(source=m.fs.pipe3.outlet,destination=m.fs.well2.inlet)

m.fs.pipe4 = pipeline(**pipeline_config,length=1600)
m.fs.pipe4.costing = idaes.core.UnitModelCostingBlock(**pipeline_costing_config)
m.fs.pipe4.diameter.fix(0.15)
m.fs.pipe4.roughness.fix(0.0475e-3)

m.fs.s8 = Arc(source=m.fs.split1.to_pipe4,destination=m.fs.pipe4.inlet)

m.fs.split3 = splitter(**splitter_config,outlet_list=['to_pipe5','to_pipe7'])

m.fs.s9 = Arc(source=m.fs.pipe4.outlet,destination=m.fs.split3.inlet)

m.fs.pipe5 = pipeline(**pipeline_config,length=800)
m.fs.pipe5.costing = idaes.core.UnitModelCostingBlock(**pipeline_costing_config)
m.fs.pipe5.diameter.fix(0.13)
m.fs.pipe5.roughness.fix(0.0475e-3)

m.fs.s10 = Arc(source=m.fs.split3.to_pipe5,destination=m.fs.pipe5.inlet)

m.fs.split4 = splitter(**splitter_config,outlet_list=['to_well3','to_pipe6'])

m.fs.s11 = Arc(source=m.fs.pipe5.outlet,destination=m.fs.split4.inlet)

m.fs.well3 = wellpad(**welltype2_config)
m.fs.well3.HCPV.fix(0.5)

m.fs.s12 = Arc(source=m.fs.split4.to_well3,destination=m.fs.well3.inlet)

m.fs.pipe6 = pipeline(**pipeline_config,length=1600)
m.fs.pipe6.costing = idaes.core.UnitModelCostingBlock(**pipeline_costing_config)
m.fs.pipe6.diameter.fix(0.1)
m.fs.pipe6.roughness.fix(0.0475e-3)

m.fs.s13 = Arc(source=m.fs.split4.to_pipe6,destination=m.fs.pipe6.inlet)

m.fs.well4 = wellpad(**welltype3_config)
m.fs.well4.HCPV.fix(0.5)

m.fs.s14 = Arc(source=m.fs.pipe6.outlet,destination=m.fs.well4.inlet)

m.fs.pipe7 = pipeline(**pipeline_config,length=800)
m.fs.pipe7.costing = idaes.core.UnitModelCostingBlock(**pipeline_costing_config)
m.fs.pipe7.diameter.fix(0.3)
m.fs.pipe7.roughness.fix(0.0475e-3)

m.fs.s15 = Arc(source=m.fs.split3.to_pipe7,destination=m.fs.pipe7.inlet)

m.fs.split5 = splitter(**splitter_config,outlet_list=['to_well5','to_pipe8'])

m.fs.s16 = Arc(source=m.fs.pipe7.outlet,destination=m.fs.split5.inlet)

m.fs.well5 = wellpad(**welltype2_config)
m.fs.well5.HCPV.fix(0.5)

m.fs.s17 = Arc(source=m.fs.split5.to_well5,destination=m.fs.well5.inlet)

m.fs.pipe8 = pipeline(**pipeline_config,length=1600)
m.fs.pipe8.costing = idaes.core.UnitModelCostingBlock(**pipeline_costing_config)
m.fs.pipe8.diameter.fix(0.13)
m.fs.pipe8.roughness.fix(0.0475e-3)

m.fs.s18 = Arc(source=m.fs.split5.to_pipe8,destination=m.fs.pipe8.inlet)

m.fs.split6 = splitter(**splitter_config,outlet_list=['to_well6','to_pipe9'])

m.fs.s19 = Arc(source=m.fs.pipe8.outlet,destination=m.fs.split6.inlet)

m.fs.well6 = wellpad(**welltype3_config)
m.fs.well6.HCPV.fix(0.5)

m.fs.s20 = Arc(source=m.fs.split6.to_well6,destination=m.fs.well6.inlet)

m.fs.pipe9 = pipeline(**pipeline_config,length=800)
m.fs.pipe9.costing = idaes.core.UnitModelCostingBlock(**pipeline_costing_config)
m.fs.pipe9.diameter.fix(0.1)
m.fs.pipe9.roughness.fix(0.0475e-3)

m.fs.s21 = Arc(source=m.fs.split6.to_pipe9,destination=m.fs.pipe9.inlet)

m.fs.well7 = wellpad(**welltype1_config)
m.fs.well7.HCPV.fix(0.5)

m.fs.s22 = Arc(source=m.fs.pipe9.outlet,destination=m.fs.well7.inlet)

#expand arcs
pyo.TransformationFactory("network.expand_arcs").apply_to(m)

#check degrees of freedom
print(f'D.o.F before specifying inlet conds={idaes.core.util.model_statistics.degrees_of_freedom(m)}')

#actual inlet specifications
m.fs.comp1.inlet.pressure[0].fix(100*100000)
m.fs.inlet_temp_constraint = pyo.Constraint(expr=m.fs.comp1.control_volume.properties_in[0].temperature==298)

#pre initialization dof
m.fs.comp1.outlet.pressure[0].fix(300*100000)
m.fs.comp1.inlet.flow_mass[0].fix(1)

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

for o in heauristic_tear_set:
    print(o.name)

def function(unit):
    try:
        unit.initialize()
    except ValueError:
        unit.deactivate_feasibility_problem()
        print(f'solve failed, propagating state')
        propagate_state(unit.inlet,unit.outlet)

seq.run(m,function)
print(f'initialization done')

#undo pre-initialization dof
m.fs.comp1.outlet.pressure[0].unfix()
m.fs.comp1.inlet.flow_mass[0].unfix()

#DoF
m.fs.comp1.inlet.flow_mass[0].unfix()
m.fs.comp1.outlet.pressure[0].unfix()

#slack variables
m.fs.pipe1.deactivate_slack_variables()
m.fs.pipe2.deactivate_slack_variables()
m.fs.pipe3.deactivate_slack_variables()
m.fs.pipe4.deactivate_slack_variables()
m.fs.pipe5.deactivate_slack_variables()
m.fs.pipe6.deactivate_slack_variables()
m.fs.pipe7.deactivate_slack_variables()
m.fs.pipe8.deactivate_slack_variables()
m.fs.pipe9.deactivate_slack_variables()

m.fs.well1.deactivate_slack_variables()
m.fs.well2.deactivate_slack_variables()
m.fs.well3.deactivate_slack_variables()
m.fs.well4.deactivate_slack_variables()
m.fs.well5.deactivate_slack_variables()
m.fs.well6.deactivate_slack_variables()
m.fs.well7.deactivate_slack_variables()

m.fs.pipe1.diameter.unfix()
m.fs.pipe2.diameter.unfix()
m.fs.pipe3.diameter.unfix()
m.fs.pipe4.diameter.unfix()
m.fs.pipe5.diameter.unfix()
m.fs.pipe6.diameter.unfix()
m.fs.pipe7.diameter.unfix()
m.fs.pipe8.diameter.unfix()
m.fs.pipe9.diameter.unfix()

#new objective function
m.fs.opex = pyo.Expression(
        expr= 365*24*3600* 4E-8*(m.fs.comp1.work_mechanical[0])
        )
m.fs.raw_mats = pyo.Expression(
    expr=365*24*3600* 0.045*m.fs.pipe1.inlet.flow_mass[0]
)
m.fs.revenue = pyo.Expression(
    expr=365*24*3600* 70*(m.fs.well1.q_OIL_PROD+m.fs.well2.q_OIL_PROD+m.fs.well3.q_OIL_PROD+m.fs.well4.q_OIL_PROD+m.fs.well5.q_OIL_PROD+m.fs.well6.q_OIL_PROD+m.fs.well7.q_OIL_PROD)
)
m.fs.capex = pyo.Expression(
    expr= m.fs.comp1.costing.capital_cost+m.fs.pipe1.costing.capital_cost+m.fs.pipe2.costing.capital_cost+m.fs.pipe3.costing.capital_cost+m.fs.pipe4.costing.capital_cost+m.fs.pipe5.costing.capital_cost+m.fs.pipe6.costing.capital_cost+m.fs.pipe7.costing.capital_cost+m.fs.pipe8.costing.capital_cost+m.fs.pipe9.costing.capital_cost
)
m.fs.obj = pyo.Objective(
    expr=(0.1*m.fs.capex+m.fs.opex+m.fs.raw_mats-m.fs.revenue)/1000000
)

m.fs.prod_target = pyo.Param(mutable=True,initialize=100000)
#constraint total production rate
m.fs.total_oil_prod_constraint = pyo.Constraint(
    expr=m.fs.prod_target/365/24/3600 == m.fs.well1.q_OIL_PROD+m.fs.well2.q_OIL_PROD+m.fs.well3.q_OIL_PROD+m.fs.well4.q_OIL_PROD+m.fs.well5.q_OIL_PROD+m.fs.well6.q_OIL_PROD+m.fs.well7.q_OIL_PROD
)

from co2_eor.util_funcs import ipopt, conopt, export_flowsheet_to_excel
import numpy as np
prod_targets = np.linspace(300000,1300000,10)

import sys
import os
from contextlib import contextmanager

# Define the suppressor
@contextmanager
def suppress_stdout():
    with open(os.devnull, "w") as devnull:
        old_stdout = sys.stdout
        sys.stdout = devnull
        try:  
            yield
        finally:
            sys.stdout = old_stdout

for prod_target in prod_targets:
    m.fs.prod_target=prod_target
    scaled_m = pyo.TransformationFactory("core.scale_model").create_using(m)
    with suppress_stdout():
        res=conopt.solve(scaled_m,tee=True)
    pyo.TransformationFactory("core.scale_model").propagate_solution(scaled_m,m)
    pyo.assert_optimal_termination(res)
    print(f'target {prod_target} solved')
    export_flowsheet_to_excel(m.fs, rf'temps/sweep_data/fs_5_{int(prod_target/1000)}.xlsx')


