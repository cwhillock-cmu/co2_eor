import pyomo.environ as pyo
import pyomo.util as pyoutil
import idaes.core
import idaes.models.properties.general_helmholtz as idaesHelmholtz
import idaes.models.unit_models.pressure_changer as idaesPressureChanger
import idaes.models.properties.swco2 as idaesSWCO2

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
            /(b.inlet_state.pressure**2-b.outlet_state.pressure**2)
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
    b.isothermal = pyo.Constraint(
            expr=b.inlet_state.temperature==b.outlet_state.temperature
            )
    #mol bal
    b.mass_bal = pyo.Constraint(
            expr=b.inlet_state.flow_mass==b.outlet_state.flow_mass
            )
    #hydraulic equation
    b.hydraulic = pyo.Constraint(
            expr=(b.inlet_state.pressure**2-b.outlet_state.pressure**2)/(2*b.average_flux_mass**2)
            ==b.average_Z*b.R/b.average_state.mw*b.average_state.temperature
            *(b.average_friction_factor*b.length/2/b.diameter+pyo.log(b.inlet_state.pressure/b.outlet_state.pressure))
            ) #isothermal compressible flow hydraulic model

    #useful expressions
    b.deltaP = pyo.Expression(
            expr=b.inlet_state.pressure-b.outlet_state.pressure
            )
    

#test block
model = pyo.ConcreteModel()
#make helmholtz parameter block
model.paramBlock =idaesHelmholtz.HelmholtzParameterBlock(
        pure_component="CO2",amount_basis=idaesHelmholtz.AmountBasis.MASS,
        state_vars=idaesHelmholtz.StateVars.TPX,phase_presentation=idaesHelmholtz.PhaseType.G,
        #has_phase_equilibrium=True, <-this does not work for some reason
        ) 
#make pipe
model.pipe1 = pyo.Block(rule = lambda b: pipe_rule(b,diameter=0.8,length=20000,roughness=0.0475e-3,paramBlock=model.paramBlock))

#fix degrees of freedom
model.pipe1.inlet_state.pressure.fix(90*100000)
model.pipe1.inlet_state.temperature.fix(293.15)
model.pipe1.inlet_state.flow_mass.fix(1266) #kg/s

#dummy objective
model.obj=pyo.Objective(expr=model.pipe1.inlet_state.pressure-model.pipe1.outlet_state.pressure)

#connect pipe to compressor
#make a new parameter block for the compressor just to get it working
model.swco2 = idaesSWCO2.SWCO2ParameterBlock()
model.compressor = idaesPressureChanger.PressureChanger(
        dynamic=False,
        property_package=model.swco2,
        compressor=True,
        thermodynamic_assumption=idaesPressureChanger.ThermodynamicAssumption.isentropic,
        )
model.compressor.display()


"""
#manually scale
scale_mult = 0.1
model.scaling_factor=pyo.Suffix(direction=pyo.Suffix.EXPORT)
model.scaling_factor[model.pipe1.diameter]=10*scale_mult
model.scaling_factor[model.pipe1.length]=1e-4*scale_mult
model.scaling_factor[model.pipe1.roughness]=1e2*scale_mult
model.scaling_factor[model.pipe1.area]=10*scale_mult
model.scaling_factor[model.pipe1.R]=1*scale_mult
model.scaling_factor[model.pipe1.inlet_velocity]=1*scale_mult
model.scaling_factor[model.pipe1.outlet_velocity]=1*scale_mult
model.scaling_factor[model.pipe1.average_velocity]=1*scale_mult
model.scaling_factor[model.pipe1.inlet_state.flow_mass]=1e-2*scale_mult
model.scaling_factor[model.pipe1.outlet_state.flow_mass]=1e-2*scale_mult
model.scaling_factor[model.pipe1.average_state.flow_mass]=1e-2*scale_mult
model.scaling_factor[model.pipe1.inlet_state.dens_mass]=1e-2*scale_mult
model.scaling_factor[model.pipe1.outlet_state.dens_mass]=1e-2*scale_mult
model.scaling_factor[model.pipe1.average_state.dens_mass]=1e-2*scale_mult
model.scaling_factor[model.pipe1.inlet_state.pressure]=1e-7*scale_mult
model.scaling_factor[model.pipe1.outlet_state.pressure]=1e-7*scale_mult
model.scaling_factor[model.pipe1.average_state.pressure]=1e-7*scale_mult
model.scaling_factor[model.pipe1.inlet_state.temperature]=1e-2*scale_mult
model.scaling_factor[model.pipe1.outlet_state.temperature]=1e-2*scale_mult
model.scaling_factor[model.pipe1.average_state.temperature]=1e-2*scale_mult
model.scaling_factor[model.pipe1.average_state.mw]=1e3*scale_mult
model.scaling_factor[model.pipe1.average_Z]=10*scale_mult
model.scaling_factor[model.pipe1.average_flux_mass]=1e-3*scale_mult
model.scaling_factor[model.pipe1.average_state.visc_d_phase["Vap"]]=1e6*scale_mult
model.scaling_factor[model.pipe1.average_reynolds]=1e-7*scale_mult
model.scaling_factor[model.pipe1.average_friction_factor]=100*scale_mult
model.scaling_factor[model.pipe1.isothermal]=1e-2*scale_mult
model.scaling_factor[model.pipe1.mass_bal]=1e4*scale_mult
model.scaling_factor[model.pipe1.hydraulic]=1e-3*scale_mult

scaled_model=pyo.TransformationFactory('core.scale_model').create_using(model)

solver = pyo.SolverFactory("ipopt")
solver.options['linear_solver']='ma57'
#solver.options['tol']=1e-6
result = solver.solve(scaled_model,logfile='ipopt_output.log',tee=True)

pyo.TransformationFactory('core.scale_model').propagate_solution(scaled_model,model)

result2=solver.solve(model,tee=True)

model.display()
print(f'pipe inlet velocity={pyo.value(model.pipe1.inlet_velocity)}')
print(f'pipe average mass density={pyo.value(model.pipe1.average_state.dens_mass)}')
print(f'pipe average Z={pyo.value(model.pipe1.average_Z)}')
print(f'mw={pyo.value(model.pipe1.average_state.mw)}')
print(f'pipe average velocity={pyo.value(model.pipe1.average_velocity)}')
print(f'pipe average reynolds={pyo.value(model.pipe1.average_reynolds)}')
print(f'pipe friction factor={pyo.value(model.pipe1.average_friction_factor)}')
print(f'pipe deltaP={pyo.value(model.pipe1.deltaP)}')
"""

