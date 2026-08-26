import pyomo.environ as pyo
import idaes.core
import pyomo.util as pyoutil
import contextlib

from idaes.models.properties.modular_properties.base.generic_property import GenericParameterBlock
from co2_eor.MPF import thermo_config
from co2_eor.gas_processing_flowsheets.aggregate_cascade import cascade

m = pyo.ConcreteModel()
m.fs = idaes.core.FlowsheetBlock(dynamic=False)
m.fs.props_vap = GenericParameterBlock(**thermo_config.configuration_vap_ideal)
m.fs.props_liq = GenericParameterBlock(**thermo_config.configuration_liq_absorption)

m.fs.column = cascade(
    vap_property_package=m.fs.props_vap,
    vap_has_phase_equilibrium=False,
    liq_property_package=m.fs.props_liq,
    liq_has_phase_equilibrium=False,
    solutes=['co2','ch4'],
    henry_coefficients={
        ('co2','A'):13.828+6.9,('co2','B'):-1720,
        ('ch4','A'):16.531+6.9,('ch4','B'):-1720,
    },
    solvents=['selexol'],
)

with open('temps/test_cascade_pprint.txt','w') as f:
    with contextlib.redirect_stdout(f):
        m.pprint()

#fix degrees of freedom
m.fs.column.top_inlet.pressure[0].fix(5*100000)
m.fs.column.top_inlet.temperature[0].fix(310)
m.fs.column.top_inlet.flow_mol_comp[0,'co2'].fix(0)
m.fs.column.top_inlet.flow_mol_comp[0,'ch4'].fix(0)
m.fs.column.top_inlet.flow_mol_comp[0,'selexol'].fix(5000)

m.fs.column.bottom_inlet.pressure[0].fix(2*100000)
m.fs.column.bottom_inlet.temperature[0].fix(300)
m.fs.column.bottom_inlet.flow_mol_comp[0,'co2'].fix(300)
m.fs.column.bottom_inlet.flow_mol_comp[0,'ch4'].fix(50)

m.fs.column.vapor_outlet_dew_point.deactivate()

print(f'DoF={idaes.core.util.model_statistics.degrees_of_freedom(m)}')

#input('paused')

m.fs.column.num_trays.fix(5)

m.fs.column.control_volume.top_in.initialize()
m.fs.column.control_volume.bot_in.initialize()
m.fs.column.control_volume.top_out.initialize()
m.fs.column.control_volume.bot_out.initialize()

print(pyo.value(m.fs.column.control_volume.top_in[0].enth_mol_phase['Liq']))
print(pyo.value(m.fs.column.control_volume.bot_in[0].enth_mol_phase['Vap']))

input('paused')

from co2_eor.util_funcs import conopt
flowsheet_solver = pyo.SolverFactory("ipopt")
flowsheet_solver.options['linear_solver']='ma97'
scaled_m = pyo.TransformationFactory("core.scale_model").create_using(m)
#solve flowsheet
res=flowsheet_solver.solve(scaled_m,tee=True)
#unscale model
pyo.TransformationFactory("core.scale_model").propagate_solution(scaled_m,m)

with open('temps/test_cascade_display.txt','w') as f:
    with contextlib.redirect_stdout(f):
        m.display()