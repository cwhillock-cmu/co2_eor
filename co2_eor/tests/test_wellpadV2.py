import pyomo.environ as pyo
import idaes.core as idaescore
import pyomo.util as pyoutil
import idaes.models.properties.general_helmholtz as idaesHelmholtz
from co2_eor import wellpad

#test block
m = pyo.ConcreteModel()
m.fs = idaescore.FlowsheetBlock(dynamic=False)
m.fs.props = idaesHelmholtz.HelmholtzParameterBlock(
        pure_component="CO2",amount_basis=idaesHelmholtz.AmountBasis.MASS,
        state_vars=idaesHelmholtz.StateVars.TPX,phase_presentation=idaesHelmholtz.PhaseType.L
        )
m.fs.wellpad = wellpad(
        property_package=m.fs.props,
        Boi=1.22,
        GOR=1.51,
        PC_A=1.386,
        PC_B=2.507,
        GB_A=0.710376,
        GB_B=1.8315,
        kovr=6.9e-12,
        SC_B=789,
        IR_base=0.001968,
        use_wellbore_pipe_model=True,
        depth=1000,
        wellbore_diameter=0.076,
        wellbore_roughness=0.0475e-3,
        )

#fix degrees of freedom
m.fs.wellpad.reservoir_state.pressure.fix(138*100000)
#m.fs.wellpad.inlet.flow_mass.fix(5)
m.fs.wellpad.reservoir_state.temperature.fix(305)
#m.fs.wellpad.control_volume.injection_state.pressure.fix(145*100000)
m.fs.wellpad.inlet.pressure.fix(74*100000)
m.fs.wellpad.inlet.temperature.fix(298)
m.fs.wellpad.HCPV.fix(0.01)

assert idaescore.util.model_statistics.degrees_of_freedom(m)==0

m.fs.wellpad.initialize()
#m.fs.wellpad.pprint()

#m.fs.wellpad.pprint()
m.obj = pyo.Objective(expr=0)
flowsheet_solver = pyo.SolverFactory("ipopt")
flowsheet_solver.options['nlp_scaling_method']='user-scaling'
res = flowsheet_solver.solve(m,tee=True)
m.display()
m.fs.wellpad.print_all()
#print(f'correction factor ={pyo.value(m.fs.wellpad.correction_factor)}')
#print(f'Oil production (STB/day)={pyo.value(m.fs.wellpad.q_OIL_PROD)*3600*24}')

