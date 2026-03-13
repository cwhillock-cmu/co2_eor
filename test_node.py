import pyomo.environ as pyo
import idaes.core as idaescore
import pyomo.util as pyoutil
import idaes.models.properties.general_helmholtz as idaesHelmholtz
import node

#test block
m = pyo.ConcreteModel()
m.fs = idaescore.FlowsheetBlock(dynamic=False)
m.fs.props = idaesHelmholtz.HelmholtzParameterBlock(
        pure_component="CO2",amount_basis=idaesHelmholtz.AmountBasis.MASS,
        state_vars=idaesHelmholtz.StateVars.TPX,phase_presentation=idaesHelmholtz.PhaseType.G,
        #has_phase_equilibrium=False,
        )

m.fs.node = node.node(
        property_package=m.fs.props,
        inlet_list=["source1","source2"],
        outlet_list=["destination1","destination2","destination3"],
        )

print(idaescore.util.model_statistics.degrees_of_freedom(m.fs.node))

#fix source1 conditions
m.fs.node.mix.source1.flow_mass[0].fix(500)
m.fs.node.mix.source1.pressure[0].fix(90*100000)
m.fs.node.mix.source1.temperature[0].fix(298)

#fix source2 conditions - cannot fix pressure because node enforces inlet pressure equality
m.fs.node.mix.source2.flow_mass[0].fix(400)
m.fs.node.mix.source2.temperature[0].fix(305)

#fix split fractions, only 2 required
m.fs.node.split.split_fraction[0,"destination1"].fix(0.6)
m.fs.node.split.split_fraction[0,"destination2"].fix(0.25)

m.fs.node.initialize()

print(idaescore.util.model_statistics.degrees_of_freedom(m.fs.node))
print(idaescore.util.model_statistics.degrees_of_freedom(m.fs.node.mix))
print(idaescore.util.model_statistics.degrees_of_freedom(m.fs.node.split))

#m.fs.node.pprint()

#print(idaescore.util.model_statistics.degrees_of_freedom(m))

solver = pyo.SolverFactory('ipopt')
res = solver.solve(m,tee=True)
m.display()
