import pyomo.environ as pyo
import idaes.core as idaescore
import pyomo.util as pyoutil
import numpy as np
from idaes.models.properties.modular_properties.base.generic_property import GenericParameterBlock
from co2_eor.thermo_config import configuration

m = pyo.ConcreteModel()
m.fs = idaescore.FlowsheetBlock(dynamic=False)
m.fs.paramBlock = GenericParameterBlock(**configuration)

m.fs.stateBlock1 = m.fs.paramBlock.build_state_block(defined_state=True)
m.pprint()
m.fs.stateBlock1.flow_mass.fix(100)
m.fs.stateBlock1.temperature.fix(298.15)

#m.fs.stateBlock1.pressure.fix(73*100000)

#m.fs.stateBlock1.initialize()
#print(f'mass density: {pyo.value(m.fs.stateBlock1.dens_mass)}')
#print(f'cp mass: {pyo.value(m.fs.stateBlock1.cp_mass)}')
solver = pyo.SolverFactory('ipopt')

plist = np.linspace(75,90,30)

for p in plist:
    m.fs.stateBlock1.pressure.fix(p*100000)
    m.fs.stateBlock1.initialize()
    res = solver.solve(m)
    print(f'pressure:{p:.2f},density:{pyo.value(m.fs.stateBlock1.dens_mass):.2f}')
    #print(f'pressure:{p:.2f},Cp: {pyo.value(m.fs.stateBlock1.cp_mol):.2f} Cp/Cv: {pyo.value(m.fs.stateBlock1.heat_capacity_ratio):.2f}')
    print(pyo.value(m.fs.stateBlock1.mw))
    print(pyo.value(m.fs.stateBlock1.enth_mol))
    print(pyo.value(m.fs.stateBlock1.visc_d_phase["Liq"]))
