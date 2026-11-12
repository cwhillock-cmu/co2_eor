import pyomo.environ as pyo
from pyomo.environ import units as pyunits
import idaes.core
from pyomo.network import Arc
import contextlib
import pickle
import time

import idaes.models.properties.general_helmholtz as idaes_helmholtz
import idaes.models.unit_models.pressure_changer as idaes_pressure_changer

from idaes.core.util.scaling import set_scaling_factor

from co2_eor import mixer, splitter, heater
from idaes.models.properties.modular_properties.base.generic_property import GenericParameterBlock
from co2_eor.liq_pipe import liqPipe
from co2_eor.wellpattern_unit_HH import wellpattern_HH
import co2_eor.MPF.thermo_config as thermo_config
from co2_eor.gas_pipe import gasPipe

from co2_eor.surrogates.rectisol_surrogate import rectisolSurrogate
from co2_eor.mpf_to_helmholtz import mpf_helmholtz_converter

from co2_eor.SSLW import SSLWCosting, SSLWCostingData

from co2_eor.utility.solvers import ipopt

output_dir = 'temps/'
m_name = 'flowsheet7_rev2_surrogate_HH'
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
    "property_package":m.fs.props_helmholtz,
    "dynamic":False,
    "compressor":True,
    "thermodynamic_assumption":idaes_pressure_changer.ThermodynamicAssumption.isentropic
}

liqpipe_config = {
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
        "momentum_mixing_type":4 #for inequality mixing type
}

splitter_config = {
        "property_package":m.fs.props_helmholtz,
        "ideal_separation":False,
        "energy_split_basis":2, #for equal_molar_enthalpy
        "momentum_balance_type":0, #for none
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

#build units
m.fs.pipe0 = liqPipe(**liqpipe_config,length=15000)
m.fs.pipe0.diameter = 1
m.fs.pipe0.roughness.fix(0.0475e-3)
m.fs.pipe0.costing = idaes.core.UnitModelCostingBlock(**liqpipe_costing_config)

m.fs.mix0 = mixer(**mixer_config,property_package=m.fs.props_helmholtz,inlet_list=['from_pipe0','from_recycle_comp'])

m.fs.main_comp = idaes_pressure_changer.PressureChanger(**compressor_config)
m.fs.main_comp.efficiency_isentropic.fix(0.85)
m.fs.main_comp.costing = idaes.core.UnitModelCostingBlock(**compressor_costing_config)

m.fs.chiller = heater(**heater_config)

m.fs.pipe1 = liqPipe(**liqpipe_config,length=2500)
m.fs.pipe1.diameter = 1
m.fs.pipe1.roughness.fix(0.0475e-3)
m.fs.pipe1.costing = idaes.core.UnitModelCostingBlock(**liqpipe_costing_config)

m.fs.split1=splitter(**splitter_config,outlet_list=['to_well1','to_pipe2'])

m.fs.well1 = wellpattern_HH(**wellpattern_config)
m.fs.well1.HCPV.fix(0.8)

m.fs.pipe2 = liqPipe(**liqpipe_config,length=3000)
m.fs.pipe2.diameter = 1
m.fs.pipe2.roughness.fix(0.0475e-3)
m.fs.pipe2.costing = idaes.core.UnitModelCostingBlock(**liqpipe_costing_config)

m.fs.well2 = wellpattern_HH(**wellpattern_config)
m.fs.well2.HCPV.fix(0.7)

m.fs.pipe3 = gasPipe(**gaspipe_config,length=3000)
m.fs.pipe3.diameter = 1
m.fs.pipe3.roughness.fix(0.0475e-3)
m.fs.pipe3.costing = idaes.core.UnitModelCostingBlock(**gaspipe_costing_config)

m.fs.mix1 = mixer(**mixer_config,property_package=m.fs.props_mix_vap_cubic,inlet_list=['from_well1','from_pipe3'])

m.fs.pipe4 = gasPipe(**gaspipe_config,length=2600)
m.fs.pipe4.diameter = 1
m.fs.pipe4.roughness.fix(0.0475e-3)
m.fs.pipe4.costing = idaes.core.UnitModelCostingBlock(**gaspipe_costing_config)

m.fs.processing_facility = rectisolSurrogate(**processingfacility_config)

m.fs.eos_converter = mpf_helmholtz_converter(
    property_package_in = m.fs.props_mix_liq_cubic,
    property_package_out = m.fs.props_helmholtz,
    conversion_type = 'convert_co2'
)

m.fs.recycle_comp = idaes_pressure_changer.PressureChanger(**compressor_config)
m.fs.recycle_comp.efficiency_isentropic.fix(0.85)
m.fs.recycle_comp.costing = idaes.core.UnitModelCostingBlock(**compressor_costing_config)

m.fs.purge_splitter = splitter(**splitter_config,outlet_list=['to_purge','to_recycle_comp'])
m.fs.purge_splitter.to_purge.flow_mol[0].lb = 0

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
m.fs.s_processing_facility_eos_converter = Arc(source=m.fs.processing_facility.recycle,destination=m.fs.eos_converter.inlet)
m.fs.s_eos_converter_purge_splitter = Arc(source=m.fs.eos_converter.outlet,destination=m.fs.purge_splitter.inlet)
m.fs.s_purge_splitter_recycle_comp = Arc(source=m.fs.purge_splitter.to_recycle_comp,destination=m.fs.recycle_comp.inlet)
m.fs.s_recycle_comp_mix0 = Arc(source=m.fs.recycle_comp.outlet,destination=m.fs.mix0.from_recycle_comp)

pyo.TransformationFactory("network.expand_arcs").apply_to(m)

#fix network specifications
m.fs.pipe0.inlet.pressure[0].fix(100*100000)
m.fs.pipe0.inlet.enth_mol[0].fix(m.fs.props_helmholtz.htpx(T=300*pyunits.K,p=100*100000*pyunits.Pa,amount_basis=idaes_helmholtz.AmountBasis.MOLE))

m.fs.chiller.outlet_temp_eq = pyo.Constraint(expr=m.fs.chiller.control_volume.properties_out[0].temperature==318)

with open(output_dir+m_name+'_postbuild_pprint.txt','w') as f:
    with contextlib.redirect_stdout(f):
        m.pprint()

#check degrees of freedom
DoF = idaes.core.util.model_statistics.degrees_of_freedom(m)
print(f'degrees of freedom={DoF}')

def initialization1():
    #set initialization dof
    m.fs.main_comp.outlet.pressure[0].fix(650*100000)
    m.fs.recycle_comp.outlet.pressure[0].fix(86*100000)
    m.fs.purge_splitter.split_fraction[0,'to_purge'].fix(1e-14)
    m.fs.purge_splitter.to_purge.pressure[0].fix(70*100000)
    #m.fs.purge_splitter.to_recycle_comp.pressure[0].fix(70*100000)
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
    m.fs.pipe0.max_pressure.fix(650*100000)
    m.fs.pipe1.max_pressure.fix(650*100000)
    m.fs.pipe2.max_pressure.fix(650*100000)
    m.fs.pipe3.max_pressure.fix(35*100000)
    m.fs.pipe4.max_pressure.fix(35*100000)

    #splitter outlet pressures
    m.fs.split1.to_well1.pressure[0].fix(600*100000)
    m.fs.split1.to_pipe2.pressure[0].fix(600*100000)

    #mixer outlet pressures
    m.fs.mix0.outlet.pressure[0].fix(84*100000)
    m.fs.mix1.outlet.pressure[0].fix(13.5*100000)

    #processing facility d.o.f
    m.fs.processing_facility.loop_solvent_flow.fix(2000)
    m.fs.processing_facility.exchanger1.outlet_temperature.fix(240)
    m.fs.processing_facility.flash3.pressure.fix(2*100000)
    m.fs.processing_facility.column.Nstages.fix(10)
    m.fs.processing_facility.deactivate_surrogate_input_bounds()

    #guess inlet flowrate
    m.fs.pipe0.inlet.flow_mol[0].fix(160)

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
        "flow_mol":{0:600},
        "enth_mol":{0:m.fs.props_helmholtz.htpx(T=300*pyunits.K,p=84*100000*pyunits.Pa,amount_basis=idaes_helmholtz.AmountBasis.MOLE)},
        "pressure":{0:84*100000},
    }
    seq.set_guesses_for(m.fs.main_comp.inlet,tear_guesses)

    def custom_initializer_func(unit):
        unit.display()
        if isinstance(unit,(liqPipe,gasPipe)):
            unit.initialize(tee=False,display_after=True)
            return
        elif isinstance(unit,wellpattern_HH):
            unit.initialize(tee=False,display_after=True)
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
    m.fs.main_comp.outlet.pressure[0].unfix()
    m.fs.recycle_comp.outlet.pressure[0].unfix()
    m.fs.purge_splitter.split_fraction[0,'to_purge'].unfix()
    m.fs.purge_splitter.to_purge.pressure[0].unfix()
    m.fs.pipe0.Pdrop_constraint.deactivate()
    m.fs.pipe1.Pdrop_constraint.deactivate()
    m.fs.pipe2.Pdrop_constraint.deactivate()
    m.fs.pipe3.Pdrop_constraint.deactivate()
    m.fs.pipe4.Pdrop_constraint.deactivate()
    m.fs.pipe0.max_pressure.unfix()
    m.fs.pipe1.max_pressure.unfix()
    m.fs.pipe2.max_pressure.unfix()
    m.fs.pipe3.max_pressure.unfix()
    m.fs.pipe4.max_pressure.unfix()

    m.fs.split1.to_well1.pressure[0].unfix()
    m.fs.split1.to_pipe2.pressure[0].unfix()

    m.fs.mix0.outlet.pressure[0].unfix()
    m.fs.mix1.outlet.pressure[0].unfix()

    #processing facility d.o.f
    m.fs.processing_facility.loop_solvent_flow.unfix()
    m.fs.processing_facility.exchanger1.outlet_temperature.unfix()
    m.fs.processing_facility.flash3.pressure.unfix()
    m.fs.processing_facility.column.Nstages.unfix()
    m.fs.processing_facility.activate_surrogate_input_bounds()

    #inlet flowrate
    m.fs.pipe0.inlet.flow_mol[0].unfix()

from co2_eor.utility.util_funcs import run_with_redirection
run_with_redirection(initialization1,output_dir+m_name+'_initialization1_log.txt')
print(f'initialization 1 done')
assert idaes.core.util.model_statistics.degrees_of_freedom(m) == DoF

#build economics
m.fs.co2_price = pyo.Param(initialize=0.068) #$/kg
m.fs.electricity_price = pyo.Param(initialize=5E-8) #$/J
m.fs.refrigeration_price = pyo.Param(initialize=3.9E-8) #$/J
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
                                    m.fs.processing_facility.costing.electrical_power
                                     +m.fs.main_comp.work_mechanical[0]
                                     +m.fs.recycle_comp.work_mechanical[0])*3600*24*365
                                  + m.fs.processing_facility.costing.refrigeration_duty*m.fs.refrigeration_price*3600*24*365)

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

with open(output_dir+m_name+'_postinitialization1_display.txt','w') as f:
    with contextlib.redirect_stdout(f):
        m.display()

from co2_eor.utility.util_funcs import report_infeasibilities
report_infeasibilities(m, output_dir+m_name+"_presolve_violations.log")

from idaes.core.scaling.autoscaling import AutoScaler
autoScaler=AutoScaler(overwrite=True)
autoScaler.scale_variables_by_magnitude(m)
#autoScaler.scale_constraints_by_jacobian_norm(m)

ipopt.options['max_iter'] = 2000
ipopt.options['linear_solver']='ma27'
ipopt.options['acceptable_tol']=1E-6
ipopt.options['tol']=1E-8
ipopt.options['expect_infeasible_problem']='yes'
#ipopt.options['nlp_scaling_method']='user-scaling'
#ipopt.options['expect_infeasible_problem']='yes'
#ipopt.options['OF_evaluate_orig_obj_at_resto_trial']='no'
#ipopt.options['OF_start_with_resto']='yes'
#scale model
scaled_m = pyo.TransformationFactory("core.scale_model").create_using(m)
res=ipopt.solve(scaled_m,tee=True)
pyo.TransformationFactory("core.scale_model").propagate_solution(scaled_m,m)

with open('temps/flowsheet7_surrogate_HH_postsolve_display.txt', 'w') as f:
    with contextlib.redirect_stdout(f):
        m.display()

report_infeasibilities(m, output_dir+m_name+"_postsolve_violations.log")

pyo.assert_optimal_termination(res)

from co2_eor.utility.util_funcs import export_flowsheet_to_excel
export_flowsheet_to_excel(m.fs, output_dir+m_name+'_surrogate_HH_solve.xlsx')

m.fs.visualize(m_name)
time.sleep(5)