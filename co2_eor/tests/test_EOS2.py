import pyomo.environ as pyo
import idaes.core as idaescore
import pyomo.util as pyoutil
import numpy as np
import contextlib
from idaes.models.properties.modular_properties.base.generic_property import GenericParameterBlock
from co2_eor.thermo_config import configuration

m = pyo.ConcreteModel()
m.fs = idaescore.FlowsheetBlock(dynamic=False)
m.fs.paramBlock = GenericParameterBlock(**configuration)
solver = pyo.SolverFactory('ipopt')

m.fs.stateBlock1 = m.fs.paramBlock.build_state_block(defined_state=False)

print(f'DoF={idaescore.util.model_statistics.degrees_of_freedom(m)}')

m.fs.stateBlock1.flow_mass.fix(100)
#m.fs.stateBlock1.flow_mol.fix(100)
m.fs.stateBlock1.temperature.fix(298.15)
m.fs.stateBlock1.pressure.fix(15*100000) #15 bar

m.fs.stateBlock1.mass_frac_comp["co2"].fix(0.8)
#m.fs.stateBlock1.mole_frac_comp["co2"].fix(0.8)

print(f'DoF={idaescore.util.model_statistics.degrees_of_freedom(m)}')

with open('temps/test_eos_stateblock_pprint_mass.txt', 'w') as f:
    with contextlib.redirect_stdout(f):
        m.pprint()

m.fs.stateBlock1.initialize()
res = solver.solve(m,tee=True)
pyo.assert_optimal_termination(res)

m.display()

print(f' mw = {pyo.value(m.fs.stateBlock1.mw)}')
print(f'dens mass = {pyo.value(m.fs.stateBlock1.dens_mass):.2f}')
print(f'enth mass = {pyo.value(m.fs.stateBlock1.enth_mass):.2f}')
print(f'visc = {pyo.value(m.fs.stateBlock1.visc_d_phase["Vap"]):.2f}')

