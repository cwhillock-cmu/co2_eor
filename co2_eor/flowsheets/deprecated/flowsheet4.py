import pyomo.environ as pyo
import pyomo.util as pyoutil
import idaes.core
import idaes.logger as idaeslog
import idaes.models.properties.general_helmholtz as idaesHelmholtz
import idaes.models.unit_models.pressure_changer as idaesPressureChanger
from co2_eor import pipeline, wellpad,node
from pyomo.network import Arc

#model, flowsheet, and parameter block
m = pyo.ConcreteModel()
m.fs = idaes.core.FlowsheetBlock(dynamic=False)
m.fs.props = idaesHelmholtz.HelmholtzParameterBlock(
        pure_component="CO2",amount_basis=idaesHelmholtz.AmountBasis.MASS,
        state_vars=idaesHelmholtz.StateVars.TPX,phase_presentation=idaesHelmholtz.PhaseType.L,
        #has_phase_equilibrium=False,
        )

#compressor 1
m.fs.comp1 = idaesPressureChanger.PressureChanger(
        property_package=m.fs.props,
        dynamic=False,
        compressor=True,
        thermodynamic_assumption=idaesPressureChanger.ThermodynamicAssumption.isentropic,
        )
m.fs.comp1.efficiency_isentropic.fix(0.85)

#pipe 1
m.fs.pipe1 = pipeline(
        property_package=m.fs.props,
        length=20000,
        diameter=0.45,
        roughness=0.0475e-3,
        alpha=5,
        ambient_temperature=293.15,
        average_pressure_type='nonlinear',
        heat_balance_type='nonisothermal',
        )

#node 1
m.fs.node1 = node(
        property_package=m.fs.props,
        inlet_list=["source1"],
        outlet_list=["branch1","branch2"]
        )

#pipe 2
m.fs.pipe2 = pipeline(
        property_package=m.fs.props,
        length=300000,
        diameter=0.4,
        roughness=0.0475e-3,
        alpha=5,
        ambient_temperature=293.15,
        average_pressure_type='nonlinear',
        heat_balance_type='nonisothermal',
        )

#pipe 3
m.fs.pipe3 = pipeline(
        property_package=m.fs.props,
        length=500,
        diameter=0.1,
        roughness=0.0473e-3,
        alpha=5,
        ambient_temperature=293.15,
        average_pressure_type='nonlinear',
        heat_balance_type='nonisothermal',
        )

#compressor 2
m.fs.comp2 = idaesPressureChanger.PressureChanger(
        property_package=m.fs.props,
        dynamic=False,
        compressor=True,
        thermodynamic_assumption=idaesPressureChanger.ThermodynamicAssumption.isentropic,
        )
m.fs.comp2.efficiency_isentropic.fix(0.86)

#node 2
m.fs.node2 = node(
        property_package=m.fs.props,
        inlet_list=["branch1"],
        outlet_list=["branch3","branch4"],
        )

#pipe 4
m.fs.pipe4 = pipeline(
        property_package=m.fs.props,
        length=600,
        diameter=0.1,
        roughness=0.0475e-3,
        alpha=5,
        ambient_temperature=293.15,
        average_pressure_type='nonlinear',
        heat_balance_type='nonisothermal',
        )

#pipe 5
m.fs.pipe5 = pipeline(
        property_package=m.fs.props,
        length=450000,
        diameter=0.3,
        roughness=0.0475e-3,
        alpha=5,
        ambient_temperature=293.15,
        average_pressure_type='nonlinear',
        heat_balance_type='nonisothermal',
        )

#wellpad 1
m.fs.wellpad1 = wellpad(
        property_package=m.fs.props,
        Boi=1.2,
        GOR=0.6,
        PC_A=0.4,
        PC_B=0.6,
        GB_A=0.62,
        GB_B=1.12,
        kovr=1.5E-14,
        use_correction_factor=False,
        BHP_max=300*100000,
        multiplier=66,
        )
m.fs.wellpad1.reservoir_state.pressure[0].fix(140*100000)
m.fs.wellpad1.reservoir_state.temperature[0].fix(306)
m.fs.wellpad1.HCPV.fix(0.5)

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

#wellpad 3
m.fs.wellpad3 = wellpad(
        property_package=m.fs.props,
        Boi=1.4,
        GOR=0.75,
        PC_A=0.25,
        PC_B=0.3,
        GB_A=0.67,
        GB_B=1.21,
        kovr=5.5E-15,
        use_correction_factor=False,
        BHP_max=340*100000,
        multiplier=100,
        )
m.fs.wellpad3.reservoir_state.pressure[0].fix(120*100000)
m.fs.wellpad3.reservoir_state.temperature[0].fix(340)
m.fs.wellpad3.HCPV.fix(0.5)

#create arcs
m.fs.s1 = Arc(source=m.fs.comp1.outlet,destination=m.fs.pipe1.inlet)
m.fs.s2 = Arc(source=m.fs.pipe1.outlet,destination=m.fs.node1.mix.source1)
m.fs.s3 = Arc(source=m.fs.node1.split.branch1,destination=m.fs.pipe2.inlet)
m.fs.s4 = Arc(source=m.fs.node1.split.branch2,destination=m.fs.pipe3.inlet)
m.fs.s5 = Arc(source=m.fs.pipe3.outlet,destination=m.fs.wellpad1.inlet)
m.fs.s6 = Arc(source=m.fs.pipe2.outlet,destination=m.fs.comp2.inlet)
m.fs.s7 = Arc(source=m.fs.comp2.outlet,destination=m.fs.node2.mix.branch1)
m.fs.s8 = Arc(source=m.fs.node2.split.branch3,destination=m.fs.pipe4.inlet)
m.fs.s9 = Arc(source=m.fs.node2.split.branch4,destination=m.fs.pipe5.inlet)
m.fs.s10 = Arc(source=m.fs.pipe4.outlet,destination=m.fs.wellpad2.inlet)
m.fs.s11 = Arc(source=m.fs.pipe5.outlet,destination=m.fs.wellpad3.inlet)

#expand arcs
pyo.TransformationFactory("network.expand_arcs").apply_to(m)

#check degrees of freedom
print(f'D.o.F before specifying inlet conds={idaes.core.util.model_statistics.degrees_of_freedom(m)}')

#inlet DoF
m.fs.comp1.inlet.pressure[0].fix(100*100000)
m.fs.comp1.inlet.temperature[0].fix(298)

#comp 1 DoF
m.fs.comp1.outlet.pressure[0].fix(350*100000)

#node 1 DoF
m.fs.node1.split.branch1.pressure[0].fix(300*100000)
m.fs.node1.split.branch2.pressure[0].fix(290*100000)

#comp 2 DoF
m.fs.comp2.outlet.pressure[0].fix(320*100000)

#node 2 DoF
m.fs.node2.split.branch3.pressure[0].fix(310*100000)
m.fs.node2.split.branch4.pressure[0].fix(320*100000)

#check degrees of freedom
print(f'D.o.F after specifying PT={idaes.core.util.model_statistics.degrees_of_freedom(m)}')

#initialize

from pyomo.network import SequentialDecomposition
seq = SequentialDecomposition()
seq.options.select_tear_method = "heuristic"
seq.options.tear_method = "Wegstein"
seq.options.iterLim = 5

G = seq.create_graph(m)
heauristic_tear_set = seq.tear_set_arcs(G,method="heuristic")
order = seq.calculation_order(G)
for o in heauristic_tear_set:
    print(o.name)

def function(unit):
    unit.initialize()

seq.run(m,function)

#activate feasibility problem
m.fs.pipe1.activate_slack_variables()
m.fs.pipe2.activate_slack_variables()
m.fs.pipe3.activate_slack_variables()
m.fs.pipe4.activate_slack_variables()
m.fs.pipe5.activate_slack_variables()
m.fs.wellpad1.activate_slack_variables()
m.fs.wellpad2.activate_slack_variables()
m.fs.wellpad3.activate_slack_variables()
m.feasibility_obj = pyo.Objective(
        expr= m.fs.pipe1.feasibility_expression + m.fs.pipe2.feasibility_expression + m.fs.pipe3.feasibility_expression
        + m.fs.pipe4.feasibility_expression + m.fs.pipe5.feasibility_expression + m.fs.wellpad1.feasibility_expression
        + m.fs.wellpad2.feasibility_expression + m.fs.wellpad3.feasibility_expression
        )

#create solver
solver1 = pyo.SolverFactory('ipopt')
solver1.options['linear_solver']='ma57'
#scale model
scaled_m = pyo.TransformationFactory("core.scale_model").create_using(m)
#solve flowsheet
res=solver1.solve(scaled_m,tee=True,logfile='ipopt_output.log')
#unscale model
pyo.TransformationFactory("core.scale_model").propagate_solution(scaled_m,m)
#m.display()
m.fs.pipe1.print_all()
m.fs.pipe2.print_all()
m.fs.pipe3.print_all()
m.fs.pipe4.print_all()
m.fs.pipe5.print_all()
m.fs.comp1.report()
m.fs.comp2.report()
m.fs.wellpad1.print_all()
m.fs.wellpad2.print_all()
m.fs.wellpad3.print_all()

input("post initialization pause")

#unfix variables for optimization

#comp 1 DoF
m.fs.comp1.outlet.pressure[0].unfix()

#node 1 DoF
m.fs.node1.split.branch1.pressure[0].unfix()
m.fs.node1.split.branch2.pressure[0].unfix()

#comp 2 DoF
m.fs.comp2.outlet.pressure[0].unfix()

#node 2 DoF
m.fs.node2.split.branch3.pressure[0].unfix()
m.fs.node2.split.branch4.pressure[0].unfix()

#deactivate feasibility problem
m.fs.pipe1.deactivate_slack_variables()
m.fs.pipe2.deactivate_slack_variables()
m.fs.pipe3.deactivate_slack_variables()
m.fs.pipe4.deactivate_slack_variables()
m.fs.pipe5.deactivate_slack_variables()
m.fs.wellpad1.deactivate_slack_variables()
m.fs.wellpad2.deactivate_slack_variables()
m.fs.wellpad3.deactivate_slack_variables()
m.feasibility_obj.deactivate()

#constraint compressor 1 work
#m.comp1Constraint = pyo.Constraint(expr=m.fs.comp1.work_mechanical[0]<=3700000)

#new objective function
m.opex_obj = pyo.Objective(
        expr= 3E-8*(m.fs.comp1.work_mechanical[0]+m.fs.comp2.work_mechanical[0])+0.050*m.fs.pipe1.inlet.flow_mass[0]-75*(m.fs.wellpad1.q_OIL_PROD+m.fs.wellpad2.q_OIL_PROD+m.fs.wellpad3.q_OIL_PROD)
        )

solver2 = pyo.SolverFactory('ipopt')
solver2.options['linear_solver']='ma97'
solver2.options['tol']=1E-6

#scale model
scaled_m = pyo.TransformationFactory("core.scale_model").create_using(m)
#solve flowsheet
res=solver2.solve(scaled_m,tee=True,logfile='ipopt_output.log')
#unscale model
pyo.TransformationFactory("core.scale_model").propagate_solution(scaled_m,m)
#m.display()

m.fs.pipe1.print_all()
m.fs.pipe2.print_all()
m.fs.pipe3.print_all()
m.fs.pipe4.print_all()
m.fs.pipe5.print_all()
m.fs.comp1.report()
m.fs.comp2.report()
m.fs.wellpad1.print_all()
m.fs.wellpad2.print_all()
m.fs.wellpad3.print_all()
"""
m.fs.pipe1.display()
m.fs.comp1.display()
m.fs.pipe2.display()
m.fs.pipe3.display()
m.fs.comp2.display()
m.fs.pipe4.display()
m.fs.pipe5.display()
m.fs.wellpad1.display()
m.fs.wellpad1.print_expressions()
m.fs.wellpad2.display()
m.fs.wellpad2.print_expressions()
m.fs.wellpad3.display()
m.fs.wellpad3.print_expressions()
"""
