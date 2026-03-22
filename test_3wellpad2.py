import pyomo.environ as pyo
import pyomo.util as pyoutil
import idaes.core
import idaes.logger as idaeslog
import idaes.models.properties.general_helmholtz as idaesHelmholtz
import idaes.models.unit_models.pressure_changer as idaesPressureChanger
from pipelineV2 import pipeline
from wellpad import wellpad
from node import node
from pyomo.network import Arc
import numpy as np

#model, flowsheet, and parameter block
m = pyo.ConcreteModel()
m.fs = idaes.core.FlowsheetBlock(dynamic=False)
m.fs.props = idaesHelmholtz.HelmholtzParameterBlock(
        pure_component="CO2",amount_basis=idaesHelmholtz.AmountBasis.MASS,
        state_vars=idaesHelmholtz.StateVars.TPX,phase_presentation=idaesHelmholtz.PhaseType.G,
        #has_phase_equilibrium=False,
        )

#wellpad 2
m.fs.wellpad2 = wellpad(
        property_package=m.fs.props,
        Boi=1.3,
        GOR=0.7,
        PC_A=0.3,
        PC_B=0.5,
        GB_A=0.53,
        GB_B=1.18,
        kovr=8.7E-15,
        use_correction_factor=False,
        BHP_max=320*100000,
        multiplier=66,
        )
m.fs.wellpad2.reservoir_state.pressure[0].fix(150*100000)
m.fs.wellpad2.reservoir_state.temperature[0].fix(295)
m.fs.wellpad2.HCPV.fix(0.5)

solver = pyo.SolverFactory('ipopt')
solver.options['linear_solver']='ma97'
solver.options['nlp_scaling_method']='user-scaling'

#d.o.f 

m.fs.wellpad2.inlet.temperature[0].fix(293.15)
injection_pressures = np.linspace(m.fs.wellpad2.reservoir_state.pressure[0].value+100000,m.fs.wellpad2.BHP_max.value,10)
IR = []
PR = []

for ip in injection_pressures:
    m.fs.wellpad2.inlet.pressure[0].fix(ip)
    m.fs.wellpad2.initialize()
    res = solver.solve(m)
    pyo.assert_optimal_termination(res)
    IR.append(pyo.value(m.fs.wellpad2.inlet.flow_mass[0]))
    PR.append(pyo.value(m.fs.wellpad2.q_OIL_PROD))

for IP, IR,PR in zip(injection_pressures,IR,PR):
    print(f'injection pressure={IP/100000}, injection rate={IR}, production rate={PR}')


