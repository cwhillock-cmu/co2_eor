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
        state_vars=idaesHelmholtz.StateVars.TPX,phase_presentation=idaesHelmholtz.PhaseType.L,
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
        "average_temperature_weight":0.5,
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
        "multiplier":6,
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
        "multiplier":4,
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
        "multiplier":4,
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

m.fs.pipe1 = pipeline(**pipeline_config,length=2000)
m.fs.pipe1.costing = idaes.core.UnitModelCostingBlock(**pipeline_costing_config)
m.fs.pipe1.diameter.fix(0.2)
m.fs.pipe1.roughness.fix(0.0475e-3)

m.fs.s1 = Arc(source=m.fs.comp1.outlet,destination=m.fs.pipe1.inlet)

m.fs.split1 = splitter(**splitter_config,outlet_list=['to_pipe2','to_pipe3'])

m.fs.s2 = Arc(source=m.fs.pipe1.outlet,destination=m.fs.split1.inlet)

m.fs.pipe2 = pipeline(**pipeline_config,length=4000)
m.fs.pipe2.costing = idaes.core.UnitModelCostingBlock(**pipeline_costing_config)
m.fs.pipe2.diameter.fix(0.15)
m.fs.pipe2.roughness.fix(0.0475e-3)

m.fs.pipe3 = pipeline(**pipeline_config,length=1000)
m.fs.pipe3.costing = idaes.core.UnitModelCostingBlock(**pipeline_costing_config)
m.fs.pipe3.diameter.fix(0.1)
m.fs.pipe3.roughness.fix(0.0475e-3)

m.fs.s3 = Arc(source=m.fs.split1.to_pipe2,destination=m.fs.pipe2.inlet)
m.fs.s4 = Arc(source=m.fs.split1.to_pipe3,destination=m.fs.pipe3.inlet)

m.fs.well1 = wellpad(**welltype3_config)
m.fs.well1.HCPV.fix(0.5)

m.fs.s5 = Arc(source=m.fs.pipe3.outlet,destination=m.fs.well1.inlet)

m.fs.split2 = splitter(**splitter_config,outlet_list=['to_well2','to_pipe4'])

m.fs.s6 = Arc(source=m.fs.pipe2.outlet,destination=m.fs.split2.inlet)

m.fs.well2 = wellpad(**welltype2_config)
m.fs.well2.HCPV.fix(0.5)

m.fs.pipe4 = pipeline(**pipeline_config,length=6000)
m.fs.pipe4.costing = idaes.core.UnitModelCostingBlock(**pipeline_costing_config)
m.fs.pipe4.diameter.fix(0.1)
m.fs.pipe4.roughness.fix(0.0475e-3)

m.fs.s7 = Arc(source=m.fs.split2.to_well2,destination=m.fs.well2.inlet)
m.fs.s8 = Arc(source=m.fs.split2.to_pipe4,destination=m.fs.pipe4.inlet)

m.fs.well3 = wellpad(**welltype1_config)
m.fs.well3.HCPV.fix(0.5)

m.fs.s9 = Arc(source=m.fs.pipe4.outlet,destination=m.fs.well3.inlet)

#expand arcs
pyo.TransformationFactory("network.expand_arcs").apply_to(m)

#check degrees of freedom
print(f'D.o.F before specifying inlet conds={idaes.core.util.model_statistics.degrees_of_freedom(m)}')

#actual inlet specifications
m.fs.comp1.inlet.pressure[0].fix(100*100000)
m.fs.comp1.inlet.temperature[0].fix(298)

#pre initialization dof
m.fs.comp1.outlet.pressure[0].fix(300*100000)
m.fs.comp1.inlet.flow_mass[0].fix(0.005)

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
print(f'pre-initialization done')

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

m.fs.well1.deactivate_slack_variables()
m.fs.well2.deactivate_slack_variables()
m.fs.well3.deactivate_slack_variables()

m.fs.pipe1.diameter.unfix()
m.fs.pipe2.diameter.unfix()
m.fs.pipe3.diameter.unfix()
m.fs.pipe4.diameter.unfix()

#m.fs.pipe1.pipe_exists.fix(1)
#m.fs.pipe2.pipe_exists.fix(1)
#m.fs.pipe3.pipe_exists.fix(0)
#m.fs.pipe4.pipe_exists.fix(1)

#constraint total production rate
m.total_oil_prod_constraint = pyo.Constraint(
    expr=200000/365/24/3600 == m.fs.well1.q_OIL_PROD+m.fs.well2.q_OIL_PROD+m.fs.well3.q_OIL_PROD
)

#new objective function
m.opex = pyo.Expression(
        expr= 365*24*3600* 4E-8*(m.fs.comp1.work_mechanical[0])
        )
m.raw_mats = pyo.Expression(
    expr=365*24*3600* 0.0*m.fs.pipe1.inlet.flow_mass[0]
)
m.revenue = pyo.Expression(
    expr=365*24*3600* 80*(m.fs.well1.q_OIL_PROD+m.fs.well2.q_OIL_PROD+m.fs.well3.q_OIL_PROD)
)
m.capex = pyo.Expression(
    expr= m.fs.comp1.costing.capital_cost+m.fs.pipe1.costing.capital_cost+m.fs.pipe2.costing.capital_cost+m.fs.pipe3.costing.capital_cost+m.fs.pipe4.costing.capital_cost
)
m.obj = pyo.Objective(
    expr=(0.1*m.capex+m.opex+m.raw_mats-m.revenue)/1000000
)

#test constraints
#m.fs.pipe3.inlet.flow_mass[0].fix(0)
#m.fs.pipe3.diameter.fix(0.0)
m.total_oil_prod_constraint.deactivate()

ipopt = pyo.SolverFactory('ipopt')
ipopt.options['linear_solver']='ma27'
ipopt.options['tol']=1E-8
ipopt.options['acceptable_tol']=1E-6
ipopt.options['halt_on_ampl_error']='yes'
bonmin = pyo.SolverFactory('bonmin')
bonmin.options['linear_solver']='ma27'
bonmin.options['tol']=1E-8
bonmin.options['acceptable_tol']=1E-6
bonmin.options['halt_on_ampl_error']='yes'
bonmin.options['bonmin.algorithm']='B-BB'
snopt = pyo.SolverFactory('snopt')
snopt.options['outlev']=2
snopt.options['major_iterations_limit'] = 10000
snopt.options['major_feasibility_tolerance']=1e-8
minos = pyo.SolverFactory('minos')
minos.options['outlev']=2
conopt=pyo.SolverFactory('conopt')
conopt.options['outlev']=3
conopt.options['logfreq']=1
conopt.options['hess']=0
baron = pyo.SolverFactory('baron')
LindoGlobal = pyo.SolverFactory('lindoglobal')
LGO = pyo.SolverFactory('LGO')
knitro = pyo.SolverFactory('knitro')
knitro.options['convex']=0
knitro.options['algorithm']=6
knitro.options['gradopt']=3
knitro.options['hessopt']=4
knitro.options['hessian_no_f']=0
knitro.options['honorbnds']=1
knitro.options['eval_fcga']=0
knitro.options['derivcheck']=0
knitro.options['cg_precond']=0
knitro.options['datacheck']=1

import contextlib
with open('temps/flowsheet_6_presolve_pprint.txt', 'w') as f:
    with contextlib.redirect_stdout(f):
        m.pprint()

#scale model
scaled_m = pyo.TransformationFactory("core.scale_model").create_using(m)
#solve flowsheet
res=bonmin.solve(scaled_m,tee=True,keepfiles=True)
#unscale model
pyo.TransformationFactory("core.scale_model").propagate_solution(scaled_m,m)

with open('temps/flowsheet_6_postsolve_display.txt', 'w') as f:
    with contextlib.redirect_stdout(f):
        m.display()

pyo.assert_optimal_termination(res)

import pandas as pd

pipe_dfs = [m.fs.pipe1.export_df(),m.fs.pipe2.export_df(),m.fs.pipe3.export_df(),m.fs.pipe4.export_df()]
well_dfs = [m.fs.well1.export_df(),m.fs.well2.export_df(),m.fs.well3.export_df()]
comp_dfs = [export_compressor_df(m.fs.comp1)]

combined_pipe_df = pd.concat(pipe_dfs)
combined_well_df = pd.concat(well_dfs)
combined_comp_df = pd.concat(comp_dfs)

flowsheet_data = {
    "pipe 1 capex":pyo.value(m.fs.pipe1.costing.capital_cost),
    "pipe 2 capex":pyo.value(m.fs.pipe2.costing.capital_cost),
    "pipe 3 capex":pyo.value(m.fs.pipe3.costing.capital_cost),
    "pipe 4 capex":pyo.value(m.fs.pipe4.costing.capital_cost),
    "compressor 1 capex":pyo.value(m.fs.comp1.costing.capital_cost),
    "total capex":pyo.value(m.capex),
    "total opex":pyo.value(m.opex),
    "total raw material cost":pyo.value(m.raw_mats),
    "total revenue":pyo.value(m.revenue),
    "objective function":pyo.value(m.obj),
}
flowsheet_df = pd.DataFrame(flowsheet_data,index=[1])

with pd.ExcelWriter('temps/flowsheet6data_localsolve.xlsx') as writer:
    combined_pipe_df.to_excel(writer,sheet_name="PipeData",index=True)
    combined_well_df.to_excel(writer,sheet_name="WellData",index=True)
    combined_comp_df.to_excel(writer,sheet_name="CompData",index=True)
    flowsheet_df.to_excel(writer,sheet_name="FlowsheetData",index=True)

