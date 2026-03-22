import pyomo.environ as pyo
import idaes.core as idaescore
import pyomo.util as pyoutil
import idaes.models.properties.general_helmholtz as idaesHelmholtz
from co2_eor import pipeline

#test block
m = pyo.ConcreteModel()
m.fs = idaescore.FlowsheetBlock(dynamic=False)
m.fs.props = idaesHelmholtz.HelmholtzParameterBlock(
        pure_component="CO2",amount_basis=idaesHelmholtz.AmountBasis.MASS,
        state_vars=idaesHelmholtz.StateVars.TPX,phase_presentation=idaesHelmholtz.PhaseType.G
        )
m.fs.pipe = pipeline.pipeline(
        property_package=m.fs.props,
        length=20000,
        diameter=0.8,
        roughness=0.0475e-3,
        )

print(idaescore.util.model_statistics.degrees_of_freedom(m))

#fix degrees of freedom
m.fs.pipe.inlet.pressure[0].fix(90*100000)
m.fs.pipe.inlet.temperature[0].fix(293.15)
m.fs.pipe.control_volume.properties_in[0].velocity.fix(3)
#m.fs.pipe.inlet.flow_mass[0].fix(1.575)
m.fs.pipe.initialize(solver='ipopt',tee=True)
m.display()
m.fs.pipe.print_expressions()
m.obj = pyo.Objective(expr=m.fs.pipe.deltaP/100000)
flowsheet_solver = pyo.SolverFactory("ipopt")
flowsheet_solver.options['nlp_scaling_method']='user-scaling'

"""
inlet_temp_list=[20,25,30,35,40,45,50,55,60]
pressure_drop_list = []
for temp in inlet_temp_list:
    m.fs.pipe.inlet.temperature.fix(temp+273.15)
    flowsheet_solver.solve(m,tee=False)
    pressure_drop_list.append(pyo.value(m.fs.pipe.deltaP)/100000)
    print(pyo.value(m.fs.pipe.feasibility_objective))

print(pressure_drop_list)
"""
