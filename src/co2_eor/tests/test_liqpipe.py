import pyomo.environ as pyo
import idaes.core as idaescore
import pyomo.util as pyoutil
import contextlib
from idaes.models.properties.modular_properties.base.generic_property import GenericParameterBlock
from co2_eor.MPF.thermo_config import configuration_liq
from co2_eor import liqPipe

#test block
m = pyo.ConcreteModel()
m.fs = idaescore.FlowsheetBlock(dynamic=False)
m.fs.props = GenericParameterBlock(**configuration_liq)

m.fs.pipe = liqPipe(
        property_package=m.fs.props,
        length=20000,
        alpha=5,
        ambient_temperature=293.15,
        average_pressure_type='nonlinear',
        heat_balance_type='nonisothermal',
        average_pressure_weight=0.5,
        average_temperature_weight=0.5,
        height_change=0,
        )

print(f'DoF={idaescore.util.model_statistics.degrees_of_freedom(m)}')

with open('temps/test_liqpipe_pprint.txt', 'w') as f:
    with contextlib.redirect_stdout(f):
        m.pprint()

#fix degrees of freedom
m.fs.pipe.diameter.fix(0.1)
m.fs.pipe.roughness.fix(0.0475e-3)
m.fs.pipe.inlet.pressure[0].fix(340*100000)
m.fs.pipe.inlet.temperature[0].fix(273.15+50)
m.fs.pipe.control_volume.properties_in[0].velocity.fix(3)
#m.fs.pipe.inlet.flow_mol_comp['co2'].fix(44000)
m.fs.pipe.inlet.flow_mol_comp[0,'ch4'].fix(1)

print(f'DoF={idaescore.util.model_statistics.degrees_of_freedom(m)}')

flowsheet_solver = pyo.SolverFactory("ipopt")
flowsheet_solver.options['linear_solver']='ma97'

m.fs.pipe.initialize(solver=flowsheet_solver,tee=True,display_after=True)

#m.obj = pyo.Objective(expr=m.fs.pipe.Pdrop**2/10000)
m.obj = pyo.Objective(expr=m.fs.pipe.diameter*10)

#scale model
scaled_m = pyo.TransformationFactory("core.scale_model").create_using(m)
#solve flowsheet
res=flowsheet_solver.solve(scaled_m,tee=True,logfile='ipopt_output.log')
#unscale model
pyo.TransformationFactory("core.scale_model").propagate_solution(scaled_m,m)

with open('temps/test_wellpatternv3_gas_pipe_display.txt', 'w') as f:
    with contextlib.redirect_stdout(f):
        m.display()

from co2_eor.util_funcs import split_print_wide_df
split_print_wide_df(m.fs.pipe.export_df())