import pyomo.environ as pyo
import idaes.core as idaescore
import pyomo.util as pyoutil
import numpy as np
import contextlib
from idaes.models.properties.modular_properties.base.generic_property import GenericParameterBlock
from co2_eor.MPF.thermo_config import configuration_liq

m = pyo.ConcreteModel()
m.fs = idaescore.FlowsheetBlock(dynamic=False)
m.fs.paramBlock = GenericParameterBlock(**configuration_liq)
solver = pyo.SolverFactory('ipopt')

m.fs.stateBlock1 = m.fs.paramBlock.build_state_block()

print(f'DoF={idaescore.util.model_statistics.degrees_of_freedom(m)}')

#m.fs.stateBlock1.flow_mass_comp['co2'].fix(100)
#m.fs.stateBlock1.flow_mass_comp['ch4'].fix(1)
m.fs.stateBlock1.flow_mol_comp['co2'].fix(100)
m.fs.stateBlock1.flow_mol_comp['ch4'].fix(1)
m.fs.stateBlock1.temperature.fix(298.15)
m.fs.stateBlock1.pressure.fix(150*100000) #150 bar

print(f'DoF={idaescore.util.model_statistics.degrees_of_freedom(m)}')

with open('temps/test_eos3_stateblock_pprint_mass.txt', 'w') as f:
    with contextlib.redirect_stdout(f):
        m.pprint()

m.fs.stateBlock1.initialize()
res = solver.solve(m,tee=True)
pyo.assert_optimal_termination(res)

m.display()

print(f' mw = {pyo.value(m.fs.stateBlock1.mw)}')
print(f'dens mass = {pyo.value(m.fs.stateBlock1.dens_mass):.2f}')
print(f'enth mass = {pyo.value(m.fs.stateBlock1.enth_mass):.2f}')
#print(f'visc = {pyo.value(m.fs.stateBlock1.visc_d_phase["Liq"])}')