import pyomo.environ as pyo
import idaes.core
import contextlib
import idaes.models.properties.modular_properties.base.generic_property
import co2_eor.MPF.thermo_config as thermo_config

m = pyo.ConcreteModel()
m.fs = idaes.core.FlowsheetBlock(dynamic=False)
m.fs.props_vap = idaes.models.properties.modular_properties.base.generic_property.GenericParameterBlock(**thermo_config.configuration_vap_cubic)
m.fs.props_liq = idaes.models.properties.modular_properties.base.generic_property.GenericParameterBlock(**thermo_config.configuration_liq_cubic)

from co2_eor.surrogates.plant_surrogate import processingFacility

surrogate_input_bounds = [
    [400,200,50e5,320,2000,260,6e5,10.5],
    [100,50,10e5,290,800,220,2e5,3.5]
]
surrogate_object = idaes.core.surrogate.AlamoSurrogate.load_from_file('temps/alamo_surrogate.json')

m.fs.plant = processingFacility(
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
        "solvent_mass_dens":791,
        "recycle_outlet_temperature":350
    },
    surrogate_object=surrogate_object
)

with open('temps/alamo_surrogate_test_pprint.txt', 'w') as f:
    with contextlib.redirect_stdout(f):
        m.pprint()

DoF = idaes.core.util.model_statistics.degrees_of_freedom(m)
print(f'degrees of freedom={DoF}')

m.fs.plant.inlet.flow_mol_comp[0,'co2'].fix(156)
m.fs.plant.inlet.flow_mol_comp[0,'ch4'].fix(84)
m.fs.plant.inlet.pressure[0].fix(2460894)
m.fs.plant.inlet.temperature[0].fix(313)
m.fs.plant.loop_solvent_flow.fix(1908)
m.fs.plant.exchanger1.outlet_temperature.fix(259)
m.fs.plant.flash3.pressure.fix(426909)
m.fs.plant.column.Nstages.fix(5)

m.fs.plant.initialize()
#m.fs.plant.custom_propagate_state()

with open('temps/alamo_surrogate_test_presolve_display.txt', 'w') as f:
    with contextlib.redirect_stdout(f):
        m.display()

from co2_eor.util_funcs import ipopt
#scale model
scaled_self = pyo.TransformationFactory('core.scale_model').create_using(m)
res = ipopt.solve(scaled_self,tee=True)
#undo scaling
pyo.TransformationFactory('core.scale_model').propagate_solution(scaled_self,m)

with open('temps/alamo_surrogate_test_postsolve_display.txt', 'w') as f:
    with contextlib.redirect_stdout(f):
        m.display()