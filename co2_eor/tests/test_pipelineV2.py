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
        state_vars=idaesHelmholtz.StateVars.TPX,phase_presentation=idaesHelmholtz.PhaseType.L
        )
m.fs.pipe = pipeline(
        property_package=m.fs.props,
        length=20000,
        alpha=5,
        ambient_temperature=293.15,
        average_pressure_type='linear',
        heat_balance_type='nonisothermal',
        average_pressure_weight=0.5,
        average_temperature_weight=0.5,
        height_change=0,
        )

#m.fs.pipe.pprint()
#m.fs.pipe.inlet_supercritical.deactivate()
#m.fs.pipe.outlet_supercritical.deactivate()

#fix degrees of freedom
m.fs.pipe.diameter.fix(0.05)
m.fs.pipe.roughness.fix(0.0475e-3)
m.fs.pipe.inlet.pressure[0].fix(340*100000)
m.fs.pipe.inlet.temperature[0].fix(273.15+50)
#m.fs.pipe.control_volume.properties_in[0].velocity.fix(3)
m.fs.pipe.inlet.flow_mass[0].fix(0)

assert idaescore.util.model_statistics.degrees_of_freedom(m)==0

flowsheet_solver = pyo.SolverFactory("ipopt")
flowsheet_solver.options['linear_solver']='ma97'

m.fs.pipe.initialize(solver=flowsheet_solver,tee=True,display_after=True)

m.obj = pyo.Objective(expr=m.fs.pipe.Pdrop**2/10000)

#scale model
scaled_m = pyo.TransformationFactory("core.scale_model").create_using(m)
#solve flowsheet
res=flowsheet_solver.solve(scaled_m,tee=True,logfile='ipopt_output.log')
#unscale model
pyo.TransformationFactory("core.scale_model").propagate_solution(scaled_m,m)

m.fs.pipe.print_all()
print(m.fs.pipe.export_df().to_string())
#print(f'mach constraint body = {pyo.value(m.fs.pipe.mach_number_constraint.body)}')
input()
inlet_temp_list=[20,25,30,35,40,45,50,55,60]
pressure_drop_list = []
T_change_list=[]
for temp in inlet_temp_list:
    m.fs.pipe.inlet.temperature.fix(temp+273.15)
    #m.fs.pipe.initialize()
    scaled_m = pyo.TransformationFactory("core.scale_model").create_using(m)
    res = flowsheet_solver.solve(scaled_m,tee=False)
    pyo.TransformationFactory("core.scale_model").propagate_solution(scaled_m,m)
    pyo.assert_optimal_termination(res)
    pressure_drop_list.append(pyo.value(m.fs.pipe.Pdrop)/100000)
    T_change_list.append(-1*pyo.value(temp+273.15-m.fs.pipe.outlet.temperature[0]))

for inlet_temp,Pdrop,outlet_temp in zip(inlet_temp_list,pressure_drop_list,T_change_list):
    print(f'inlet_temp: {inlet_temp}, P drop: {Pdrop}, temp change: {outlet_temp}')

