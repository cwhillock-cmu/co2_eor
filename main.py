import pyomo.environ as pyo
import pyomo.util as pyoutil
import idaes.core
import idaes.models.properties.general_helmholtz as idaesHelmholtz
from idaes.models.properties.swco2 import htpx

#define the pipe block (in a function because thats what I know how to do rn)
def pipe_rule(b,diameter,length,roughness,paramBlock):
    #declare parameters
    b.diameter = pyo.Param(initialize=diameter)#m
    b.length = pyo.Param(initialize=length) #m
    b.roughness = pyo.Param(initialize=roughness) #m
    b.area = pyo.Expression(
            expr=3.1415926*b.diameter**2/4
            ) #m^2
    b.R = pyo.Param(initialize=8.314462) #m^3*Pa*K^{-1}*mol^{-1}
    
    #get thermodynamic state blocks
    b.inlet_state = paramBlock.build_state_block(defined_state=False)
    b.outlet_state = paramBlock.build_state_block(defined_state=False)
    b.average_state = paramBlock.build_state_block(defined_state=False)

    #inlet variables
    b.inlet_velocity = pyo.Expression(
            expr=b.inlet_state.flow_mass/b.inlet_state.dens_mass/b.area
            )

    #outlet variables
    b.outlet_velocity = pyo.Expression(
            expr=b.outlet_state.flow_mass/b.outlet_state.dens_mass/b.area
            )

    #average state variables
    b.average_pressure_constraint = pyo.Constraint(
            expr=b.average_state.pressure*3==
            2*(b.inlet_state.pressure**3-b.outlet_state.pressure**3)
            /(b.inlet_state.pressure**2-b.outlet_state.pressure**3)
            )
    b.average_temperature_constraint = pyo.Constraint(
            expr=b.average_state.temperature*2==
            b.inlet_state.temperature+b.outlet_state.temperature
            )
    b.average_flow_mass_constraint = pyo.Constraint(
            expr=b.average_state.flow_mass*2==
            b.inlet_state.flow_mass+b.outlet_state.flow_mass
            ) #only at steady state
    b.average_velocity = pyo.Expression(
            expr=b.average_state.flow_mass/b.average_state.dens_mass/b.area
            )
    b.average_Z = pyo.Expression(
            expr=b.average_state.pressure*b.average_state.mw
            /b.average_state.dens_mass/b.R/b.average_state.temperature
            )
    b.average_flux_mass = pyo.Expression(
            expr=b.average_state.dens_mass*b.average_velocity
            )
    b.average_reynolds = pyo.Expression(
            expr=b.average_flux_mass*b.diameter/b.average_state.visc_d_phase["Vap"]
            ) #lock this in to vapor phase for now, add modularity later
    b.average_friction_factor = pyo.Expression(
            expr=0.25*
            (pyo.log10(b.roughness/3.7/b.diameter+5.74/(b.average_reynolds**0.9)))**(-2)
            )#swamee-jain friction factor

    #constraints
    #isothermal
    b.isothermal_constraint = pyo.Constraint(
            expr=b.inlet_state.temperature==b.outlet_state.temperature
            )
    #mol bal
    b.mass_bal = pyo.Constraint(
            expr=b.inlet_state.flow_mass==b.outlet_state.flow_mass
            )
    #hydraulic equation
    b.hydraulic = pyo.Constraint(
            expr=(b.inlet_state.pressure**2-b.outlet_state.pressure**2)/(2*b.average_flux_mass**2)
            ==b.average_Z*b.R*b.average_state.mw*b.average_state.temperature
            *(b.average_friction_factor*b.length/2/b.diameter+pyo.log(b.inlet_state.pressure/b.outlet_state.pressure))
            ) #isothermal compressible flow hydraulic model
    

#test block
model = pyo.ConcreteModel()
#make helmholtz parameter block
model.paramBlock =idaesHelmholtz.HelmholtzParameterBlock(
        pure_component="CO2",amount_basis=idaesHelmholtz.AmountBasis.MASS,
        state_vars=idaesHelmholtz.StateVars.TPX,phase_presentation=idaesHelmholtz.PhaseType.G,
        ) 
#make pipe
model.pipe1 = pyo.Block(rule = lambda b: pipe_rule(b,diameter=0.8,length=50000,roughness=1.0475e-3,paramBlock=model.paramBlock))

#check model
print(idaes.core.util.model_statistics.degrees_of_freedom(model))

#fix degrees of freedom
model.pipe1.inlet_state.pressure.fix(118*100000*pyo.units.Pa)
model.pipe1.inlet_state.temperature.fix(301.15*pyo.units.K)
model.pipe1.inlet_state.flow_mass.fix(20*1000/365/24/3600) #kg/s

#dummy objective
model.obj = pyo.Objective(expr=1)

solver = pyo.SolverFactory("ipopt")
result = solver.solve(model,logfile='ipopt_output.log',tee=True)
print(model.display())

