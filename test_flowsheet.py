import pyomo.environ as pyo
import pyomo.util as pyoutil
import idaes.core
import idaes.logger as idaeslog
import idaes.models.properties.general_helmholtz as idaesHelmholtz
import idaes.models.unit_models.pressure_changer as idaesPressureChanger
from pipeline import pipeline
from wellpad import wellpad
from pyomo.network import Arc
 
#create flowsheet
model = pyo.ConcreteModel()
model.flowsheet = idaes.core.FlowsheetBlock(
        dynamic=False
        )
#make helmholtz parameter block
model.flowsheet.paramBlock =idaesHelmholtz.HelmholtzParameterBlock(
        pure_component="CO2",amount_basis=idaesHelmholtz.AmountBasis.MASS,
        state_vars=idaesHelmholtz.StateVars.TPX,phase_presentation=idaesHelmholtz.PhaseType.G,
        #has_phase_equilibrium=True, <-this does not work for some reason
        ) 
#make pipe
model.flowsheet.pipe1 = pipeline(
        property_package=model.flowsheet.paramBlock,
        length=20000,
        diameter = 0.8,
        roughness=0.0475e-3,
        )

#make wellpad
model.flowsheet.wellpad1 = wellpad(
        property_package=model.flowsheet.paramBlock,
        Boi=1.22,
        GOR=1.51,
        PC_A=1.386,
        PC_B=2.507,
        GB_A=0.710376,
        GB_B=1.8315,
        kovr=6.9e-13,
        SC_B=789,
        IR_base=0.001968,
        )

#create compressor
model.flowsheet.compressor = idaesPressureChanger.PressureChanger(
        dynamic=False,
        property_package=model.flowsheet.paramBlock,
        compressor=True,
        thermodynamic_assumption=idaesPressureChanger.ThermodynamicAssumption.isentropic,
        )

#fix pipe degrees of freedom
model.flowsheet.pipe1.inlet.pressure[0].fix(90*100000)
model.flowsheet.pipe1.inlet.temperature[0].fix(298)

#create constraints to connect pipe and compressor
model.flowsheet.s1 = Arc(source=model.flowsheet.pipe1.outlet,destination=model.flowsheet.compressor.inlet)

#define other compressor constraints
model.flowsheet.compressor.outlet.pressure[0].fix(140*100000)
model.flowsheet.compressor.efficiency_isentropic.fix(0.85)

#fix wellpad degrees of freedom
model.flowsheet.wellpad1.reservoir_state.pressure[0].fix(138*100000)
model.flowsheet.wellpad1.reservoir_state.temperature[0].fix(306)
model.flowsheet.wellpad1.HCPV.fix(0.3)

#connect compressor to wellpad
model.flowsheet.s2 = Arc(source=model.flowsheet.compressor.outlet,destination=model.flowsheet.wellpad1.inlet)

#expand arcs
pyo.TransformationFactory('network.expand_arcs').apply_to(model)

#initialize flowsheet and units
model.flowsheet.pipe1.inlet.flow_mass[0].fix(0.3)
model.flowsheet.compressor.outlet.pressure[0].unfix()
assert idaes.core.util.model_statistics.degrees_of_freedom(model)==0
model.flowsheet.pipe1.initialize()
model.flowsheet.compressor.initialize(outlvl=idaeslog.INFO)
model.flowsheet.wellpad1.initialize()
model.flowsheet.pipe1.inlet.flow_mass[0].unfix()
model.flowsheet.compressor.outlet.pressure[0].fix(140*100000)


#solve feasibility problem
model.flowsheet.pipe1.activate_feasibility_problem()
model.flowsheet.wellpad1.activate_feasibility_problem()
model.flowsheet.pipe1.feasibility_objective.deactivate()
model.flowsheet.wellpad1.feasibility_objective.deactivate()
model.obj = pyo.Objective(expr=model.flowsheet.pipe1.feasibility_expression + model.flowsheet.wellpad1.feasibility_expression)
solver = pyo.SolverFactory("ipopt")
solver.options['linear_solver']='ma27'
solver.options['nlp_scaling_method']='user-scaling'

res = solver.solve(model,tee=True)

model.flowsheet.pipe1.display()
model.flowsheet.compressor.display()
model.flowsheet.wellpad1.display()
model.flowsheet.wellpad1.print_expressions()

