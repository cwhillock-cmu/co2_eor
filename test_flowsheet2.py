import pyomo.environ as pyo
import pyomo.util as pyoutil
import idaes.core
import idaes.logger as idaeslog
import idaes.models.properties.general_helmholtz as idaesHelmholtz
import idaes.models.unit_models.pressure_changer as idaesPressureChanger
from pipeline import pipeline
from wellpad import wellpad
from node import node
from pyomo.network import Arc
 
#create flowsheet
m = pyo.ConcreteModel()
m.fs = idaes.core.FlowsheetBlock(
        dynamic=False
        )
#make helmholtz parameter block
m.fs.paramBlock =idaesHelmholtz.HelmholtzParameterBlock(
        pure_component="CO2",amount_basis=idaesHelmholtz.AmountBasis.MASS,
        state_vars=idaesHelmholtz.StateVars.TPX,phase_presentation=idaesHelmholtz.PhaseType.G,
        #has_phase_equilibrium=True, <-this does not work for some reason
        )

#create pipes
m.fs.pipe1 = pipeline(
        property_package=m.fs.paramBlock,
        length=20000, #m
        diameter=0.8, #m
        roughness=0.0475e-3, #m
        )
m.fs.pipe2 = pipeline(
        property_package=m.fs.paramBlock,
        length=20000, #m
        diameter= 0.8, #m
        roughness = 0.0475e-3, #m
        )
m.fs.pipe3 = pipeline(
        property_package=m.fs.paramBlock,
        length=13000, #m
        diameter = 0.8, #m
        roughness = 0.0475e-3, #m
        )
m.fs.pipe4 = pipeline(
        property_package=m.fs.paramBlock,
        length=15000, #m
        diameter = 0.8, #m
        roughness = 0.0475e-3, #m
        )
m.fs.pipe5 = pipeline(
        property_package=m.fs.paramBlock,
        length=17000, #m
        diameter=0.8, #m
        roughness = 0.0475e-3, #m
        )
#create compressor
m.fs.compressor1 = idaesPressureChanger.PressureChanger(
        property_package=m.fs.paramBlock,
        dynamic=False,
        compressor=True,
        thermodynamic_assumption=idaesPressureChanger.ThermodynamicAssumption.isentropic,
        )
#fix compressor specifications - so I don't forget to do it later
m.fs.compressor1.efficiency_isentropic.fix(0.85)
#create node
m.fs.node1 = node(
        property_package=m.fs.paramBlock,
        inlet_list=["source1","source2"],
        outlet_list=["destination1","destination2","destination3"]
        )

#connect units
m.fs.s1 = Arc(source=m.fs.pipe1.outlet,destination=m.fs.node1.mix.source1)
m.fs.s2 = Arc(source=m.fs.compressor1.outlet,destination=m.fs.pipe2.inlet)
m.fs.s3 = Arc(source=m.fs.pipe2.outlet,destination=m.fs.node1.mix.source2)
m.fs.s4 = Arc(source=m.fs.node1.split.destination1,destination=m.fs.pipe3.inlet)
m.fs.s5 = Arc(source=m.fs.node1.split.destination2,destination=m.fs.pipe4.inlet)
m.fs.s6 = Arc(source=m.fs.node1.split.destination3,destination=m.fs.pipe5.inlet)
#expand arcs
pyo.TransformationFactory("network.expand_arcs").apply_to(m)

#check degrees of freedom
print(f'D.o.F before specifying inlet/outlet conds={idaes.core.util.model_statistics.degrees_of_freedom(m)}')

#specify inlet and outlet conditions
#pipe 1
m.fs.pipe1.inlet.pressure[0].fix(200*100000)
m.fs.pipe1.inlet.temperature[0].fix(305)
#compressor 1
m.fs.compressor1.inlet.pressure[0].fix(130*100000)
m.fs.compressor1.inlet.temperature[0].fix(273)
#pipe 3
m.fs.pipe3.outlet.pressure[0].fix(150*100000)
#pipe 4
m.fs.pipe4.outlet.pressure[0].fix(145*100000)
#pipe 5
m.fs.pipe5.outlet.pressure[0].fix(140*100000)

print(f'D.o.F after specifying inlet/outlet conds={idaes.core.util.model_statistics.degrees_of_freedom(m)}')

#specify split fractions
m.fs.node1.split.split_fraction[0,"destination1"].fix(0.25)
m.fs.node1.split.split_fraction[0,"destination2"].fix(0.4)

print(f'D.o.F after specifying split fractions={idaes.core.util.model_statistics.degrees_of_freedom(m)}')

#specify inlet flowrates
m.fs.pipe1.inlet.flow_mass[0].fix(500)
m.fs.compressor1.inlet.flow_mass[0].fix(425)

print(f'D.o.F after specifying inlet flowrates={idaes.core.util.model_statistics.degrees_of_freedom(m)}')

#initialize units
m.fs.pipe1.initialize()
#input('paused')
m.fs.compressor1.ratioP.fix(1.1)
m.fs.compressor1.initialize()
m.fs.compressor1.ratioP.unfix()
#input('paused')
m.fs.pipe2.initialize()
#input('paused')
m.fs.node1.initialize()
#input('paused')
m.fs.pipe3.initialize()
#input('paused')
m.fs.pipe4.initialize()
#input('paused')
m.fs.pipe5.initialize()

#activate feasibility objective
m.fs.pipe1.activate_feasibility_problem()
m.fs.pipe1.feasibility_objective.deactivate()
m.fs.pipe2.activate_feasibility_problem()
m.fs.pipe2.feasibility_objective.deactivate()
m.fs.pipe3.activate_feasibility_problem()
m.fs.pipe3.feasibility_objective.deactivate()
m.fs.pipe4.activate_feasibility_problem()
m.fs.pipe4.feasibility_objective.deactivate()
m.fs.pipe5.activate_feasibility_problem()
m.fs.pipe5.feasibility_objective.deactivate()
m.obj = pyo.Objective(
        expr= m.fs.pipe1.feasibility_expression + m.fs.pipe2.feasibility_expression
        + m.fs.pipe3.feasibility_expression +m.fs.pipe4.feasibility_expression
        + m.fs.pipe5.feasibility_expression
        )
#scale model
scaled_m = pyo.TransformationFactory("core.scale_model").create_using(m)
#create solver
solver = pyo.SolverFactory('ipopt')
solver.options['linear_solver']='ma97'
#solve flowsheet
res=solver.solve(scaled_m,tee=True)
#unscale model
pyo.TransformationFactory("core.scale_model").propagate_solution(scaled_m,m)
m.display()


