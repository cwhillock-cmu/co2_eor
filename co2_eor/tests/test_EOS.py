import pyomo.environ as pyo
import idaes.core as idaescore
import pyomo.util as pyoutil
import idaes.models.properties.general_helmholtz as idaesHelmholtz
import numpy as np

m = pyo.ConcreteModel()
m.fs = idaescore.FlowsheetBlock(dynamic=False)
m.fs.paramBlock = idaesHelmholtz.HelmholtzParameterBlock(
        pure_component="CO2",amount_basis=idaesHelmholtz.AmountBasis.MASS,
        state_vars=idaesHelmholtz.StateVars.TPX,phase_presentation=idaesHelmholtz.PhaseType.MIX
        )

m.fs.stateBlock1 = m.fs.paramBlock.build_state_block(defined_state=True)

m.fs.stateBlock1.flow_mass.fix(100)
m.fs.stateBlock1.temperature.fix(293)
#m.fs.stateBlock1.pressure.fix(73*100000)

#m.fs.stateBlock1.initialize()
#print(f'mass density: {pyo.value(m.fs.stateBlock1.dens_mass)}')
#print(f'cp mass: {pyo.value(m.fs.stateBlock1.cp_mass)}')

plist = np.linspace(60,120,40)

for p in plist:
    m.fs.stateBlock1.pressure=p*100000
    m.fs.stateBlock1.vapor_frac=0.5
    m.fs.stateBlock1.initialize()
    print(f'pressure:{p:.2f},density:{pyo.value(m.fs.stateBlock1.dens_mass):.2f},Tcrit:{pyo.value(m.fs.stateBlock1.temperature_crit):.2f},Pcrit:{pyo.value(m.fs.stateBlock1.pressure_crit)/100000:.2f},Tsat:{pyo.value(m.fs.stateBlock1.temperature_sat):.2f},Psat:{pyo.value(m.fs.stateBlock1.pressure_sat)/100000:.2f},L Visc:{pyo.value(m.fs.stateBlock1.visc_d_phase["Liq"])},G Visc:{pyo.value(m.fs.stateBlock1.visc_d_phase["Vap"])},x:{pyo.value(m.fs.stateBlock1.vapor_frac)}')
