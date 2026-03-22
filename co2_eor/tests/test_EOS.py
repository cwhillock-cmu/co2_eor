import pyomo.environ as pyo
import idaes.core as idaescore
import pyomo.util as pyoutil
import idaes.models.properties.general_helmholtz as idaesHelmholtz
import numpy as np

m = pyo.ConcreteModel()
m.fs = idaescore.FlowsheetBlock(dynamic=False)
m.fs.paramBlock = idaesHelmholtz.HelmholtzParameterBlock(
        pure_component="CO2",amount_basis=idaesHelmholtz.AmountBasis.MASS,
        state_vars=idaesHelmholtz.StateVars.TPX,phase_presentation=idaesHelmholtz.PhaseType.G
        )

m.fs.stateBlock1 = m.fs.paramBlock.build_state_block(defined_state=True)

m.fs.stateBlock1.flow_mass=100
m.fs.stateBlock1.temperature=283.15
m.fs.stateBlock1.pressure=73*100000

m.fs.stateBlock1.initialize()
print(f'mass density: {pyo.value(m.fs.stateBlock1.dens_mass)}')
print(f'cp mass: {pyo.value(m.fs.stateBlock1.cp_mass)}')

plist = np.linspace(70,80,20)

for p in plist:
    m.fs.stateBlock1.pressure=p*100000
    m.fs.stateBlock1.initialize()
    print(f'pressure:{p},density:{pyo.value(m.fs.stateBlock1.dens_mass)}')
