import pyomo.environ as pyo
import idaes.core
import contextlib
import idaes.models.properties.modular_properties.base.generic_property
import co2_eor.MPF.thermo_config as thermo_config

m = pyo.ConcreteModel()
m.fs = idaes.core.FlowsheetBlock(dynamic=False)
m.fs.props_vap = idaes.models.properties.modular_properties.base.generic_property.GenericParameterBlock(**thermo_config.configuration_vap_cubic)
m.fs.props_liq = idaes.models.properties.modular_properties.base.generic_property.GenericParameterBlock(**thermo_config.configuration_liq_cubic)

from co2_eor.surrogates.rectisol_surrogate import rectisolSurrogate

surrogate_input_bounds = [
    [400,200,50e5,320,2000,260,6e5,10.5],
    [100,50,10e5,290,800,220,2e5,3.5]
]
hyperplane_list = [
    {"weights": [1.12241413e-03, -2.56485953e-02, 1.67543279e-06, 2.40731609e-03, 1.42386045e-03, 1.04199595e-02, -6.81482689e-06, 9.81953962e-01],
     "bias":-15.0399895}
]
surrogate_object = idaes.core.surrogate.AlamoSurrogate.load_from_file('temps/alamo_surrogate.json')

m.fs.plant = rectisolSurrogate(
    inlet_property_package=m.fs.props_vap,
    outlet_property_package=m.fs.props_vap,
    recycle_property_package=m.fs.props_liq,
    surrogate_input_bounds=surrogate_input_bounds,
    unit_model_parameters={
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
    global_parameters={
        "solvent_MW":32.04e-3,
        "solvent_mass_dens":800,
        "recycle_outlet_temperature":350
    },
    surrogate_object=surrogate_object,
    hyperplane_cuts=hyperplane_list,
)

with open('temps/alamo_surrogate_test_pprint.txt', 'w') as f:
    with contextlib.redirect_stdout(f):
        m.pprint()

DoF = idaes.core.util.model_statistics.degrees_of_freedom(m)
print(f'degrees of freedom={DoF}')

m.fs.plant.inlet.flow_mol_comp[0,'co2'].fix(200)
m.fs.plant.inlet.flow_mol_comp[0,'ch4'].fix(100)
m.fs.plant.inlet.pressure[0].fix(25e5)
m.fs.plant.inlet.temperature[0].fix(300)
m.fs.plant.loop_solvent_flow.fix(1900)
m.fs.plant.exchanger1.outlet_temperature.fix(250)
m.fs.plant.flash3.pressure.fix(2e5)
#m.fs.plant.column.Nstages.fix(5)

m.fs.plant.initialize()
#m.fs.plant.custom_propagate_state()

with open('temps/alamo_surrogate_test_presolve_display.txt', 'w') as f:
    with contextlib.redirect_stdout(f):
        m.display()

from co2_eor.utility.solvers import ipopt
#scale model
scaled_self = pyo.TransformationFactory('core.scale_model').create_using(m)
res = ipopt.solve(scaled_self,tee=True)
#undo scaling
pyo.TransformationFactory('core.scale_model').propagate_solution(scaled_self,m)

with open('temps/alamo_surrogate_test_postsolve_display.txt', 'w') as f:
    with contextlib.redirect_stdout(f):
        m.display()