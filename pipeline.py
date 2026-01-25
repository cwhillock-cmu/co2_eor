import pyomo.environ as pyo
import pyomo.util as pyoutil
import idaes.core
import idaes.models.properties.general_helmholtz as idaesHelmholtz
from idaes.models.properties.swco2 import htpx

pipe = pyo.ConcreteModel()

#declare parameters (design variables)
pipe.dia = pyo.Param(initialize=0.8) #m
pipe.length = pyo.Param(initialize=50000) #m
pipe.roughness = pyo.Param(initialize=1.0475e-3) #m
pipe.area = pyo.Param(initialize=3.1415926*pipe.dia**2/4) #m^2

#declare global variables


#get thermodynamic stuff
pipe.paramBlock = idaesHelmholtz.HelmholtzParameterBlock(
        pure_component="CO2",amount_basis=idaesHelmholtz.AmountBasis.MOLE,
        state_vars=idaesHelmholtz.StateVars.PH,phase_presentation=idaesHelmholtz.PhaseType.G,
        )
pipe.inlet_state = pipe.paramBlock.build_state_block(defined_state=False)
pipe.outlet_state = pipe.paramBlock.build_state_block(defined_state=False)

#inlet variables
pipe.inlet_pressure = pyo.Var(initialize=118*100000) #Pa
pipe.inlet_pressure.fix()

pipe.inlet_temperature = pyo.Var(initialize=28+273.15) #K
pipe.inlet_temperature.fix()

pipe.inlet_pressure_constrain_to_stateblock = pyo.Constraint(
        expr=pipe.inlet_state.pressure==pipe.inlet_pressure
        )
pipe.inlet_enthalpy_calculation_constraint = pyo.Constraint(
        expr=pipe.inlet_state.enthalpy==htpx(T=pipe.inlet_temperature*pyo.units.K,P=pipe.inlet_pressure*pyo.units.Pa)
        )

pipe.inlet_velocity= pyo.Expression(
        expr=pipe.inlet_state.dens_mol*pipe.area
        )


#outlet variables
pipe.outlet_pressure = pyo.Var() #Pa
pipe.outlet_temperature = pyo.Var() #K
pipe_outlet_state.pressure=pipe.outlet_pressure #Pa
pipe.outlet_state.enthalpy = htpx(T=pipe.outlet_temperature*pyo.units.K,P=pipe.outlet_pressure*pyo.units.Pa)

#average variables(use expression?)
pipe.average_pressure = pyo.Expression(
        expr=2/3*(pipe.inlet_pressure**3-pipe.outlet_pressure**3)/(pipe.inlet_pressure**2-pipe.outlet_pressure**2)
        )


#constraints
#isothermal
pipe.isothermal = pyo.Constraint(
        expr=pipe.inlet_temperature==pipe.outlet_temperature
        )
#mass bal
pipe.mass_bal = pyo.Constraint(
        expr=pipe.inlet_state.flow_mol==pipe.outlet_state.flow_mol
        )


"""
h = htpx(T=298.15*pyo.units.K,P=120*10000*pyo.units.Pa)
print(f'enthalpy test={h}')
pipe.stateBlock.display()
pipe.stateBlock.report()
pipe.stateBlock.flow_mol.fix(100)
pipe.stateBlock.pressure.fix(120)
pipe.stateBlock.enth_mol.fix(h)
"""

