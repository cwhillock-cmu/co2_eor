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
m.fs.pipe = pipeline(
        property_package=m.fs.props,
        length=20000,
        diameter=0.8,
        roughness=0.0475e-3,
        alpha=5,
        ambient_temperature=293.15,
        average_pressure_type='nonlinear',
        heat_balance_type='nonisothermal',
        )

print(idaescore.util.model_statistics.degrees_of_freedom(m))
#m.fs.pipe.pprint()

#fix degrees of freedom
m.fs.pipe.inlet.pressure[0].fix(90*100000)
m.fs.pipe.inlet.temperature[0].fix(313.15)
m.fs.pipe.control_volume.properties_in[0].velocity.fix(3)
#m.fs.pipe.inlet.flow_mass[0].fix(40)

flowsheet_solver = pyo.SolverFactory("bonmin")
#flowsheet_solver.options['linear_solver']='ma97'
print(flowsheet_solver.options)

m.fs.pipe.initialize(solver=flowsheet_solver,tee=True,display_after=True)

m.obj = pyo.Objective(expr=m.fs.pipe.Pdrop/100000)

#scale model
scaled_m = pyo.TransformationFactory("core.scale_model").create_using(m)
#solve flowsheet
res=flowsheet_solver.solve(scaled_m,tee=True,logfile='ipopt_output.log')
#unscale model
pyo.TransformationFactory("core.scale_model").propagate_solution(scaled_m,m)

"""
m.display()
m.fs.pipe.print_all()
print(pyo.value(m.fs.pipe.control_volume.properties_avg.cp_mass))
print(pyo.value(m.fs.pipe.control_volume.properties_avg.dens_mass))
print(pyo.value(m.fs.pipe.control_volume.properties_avg.visc_d_phase["Vap"]))
print(pyo.value(m.fs.pipe.control_volume.properties_avg.mw))
print()
print(pyo.value(m.fs.pipe.control_volume.properties_out[0].dens_mass))
print(m.fs.pipe.name)
"""
m.fs.pipe.print_all()
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
