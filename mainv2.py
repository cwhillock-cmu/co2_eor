import pyomo.environ as pyo
import pyomo.util as pyoutil
import idaes.core
import idaes.logger as idaeslog
import idaes.models.properties.general_helmholtz as idaesHelmholtz
import idaes.models.unit_models.pressure_changer as idaesPressureChanger


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

def reservoir_rule(b,reservoir_pressure,reservoir_temperature,BOI,HCPV,PC_A,PC_B,GB_A,GB_B,GOR,paramBlock):
    #no hydraulics yet

    #declare parameters
    b.BOI = pyo.Param(initialize=BOI) #RB/STB
    b.PC_A = pyo.Param(initialize=PC_A) #production curve fit parameter a
    b.PC_B = pyo.Param(initialize=PC_B) #production curve fit parameter b
    b.GB_A = pyo.Param(initialize=GB_A) #gas breakthrough curve fit parameter a
    b.GB_B = pyo.Param(initialize=GB_B) #gas breakthrough curve fit parameter b
    b.GOR = pyo.Param(initialize=GOR) #gas produced per oil produced (not from gas breakthrough) #RB/STB

    #state blocks
    b.inlet_state = paramBlock.build_state_block(defined_state=False)
    b.reservoir_state = paramBlock.build_state_block(defined_state=False)
    #fix pressure and temperature of reservoir state
    b.reservoir_state.pressure.fix(reservoir_pressure) #Pa
    b.reservoir_state.temperature.fix(reservoir_temperature) #K
    #mass balance between inlet state and reservoir state
    b.mass_bal = pyo.Constraint(
            expr=b.inlet_state.flow_mass==b.reservoir_state.flow_mass
            )

    #HCPV is the cumulative volume of CO2 injected into the well, for now it will just be a parameter
    b.HCPV = pyo.Param(initialize=HCPV) #unitless

    #CO2 injection rate
    b.V_CO2_inj = pyo.Expression(
            expr=b.reservoir_state.flow_mass/b.reservoir_state.dens_mass
            ) #RB/s
    #oil production rate
    b.V_oil_prod = pyo.Expression(
            expr=1/b.BOI*b.PC_A*b.PC_B/(b.HCPV+b.PC_B)**2*b.V_CO2_inj
            ) #STB/s
    #gas production rate
    b.V_gas_prod = pyo.Expression(
            expr=b.GOR*b.V_oil_prod*b.BOI
            ) #RB/s
    #gas breakthrough rate
    b.V_CO2_brkth = pyo.Expression(
            expr=b.GB_A*b.GB_B*b.HCPV**(b.GB_B-1)*b.V_CO2_inj
            ) #RB/s
    b.V_gas_out = pyo.Expression(
            expr=b.V_gas_prod + b.V_CO2_brkth
            )

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
model.flowsheet.pipe1 = pyo.Block(rule = lambda b: pipe_rule(b,diameter=0.8,length=20000,roughness=0.0475e-3,paramBlock=model.flowsheet.paramBlock))

#fix degrees of freedom
model.flowsheet.pipe1.inlet_state.pressure.fix(90*100000)
model.flowsheet.pipe1.inlet_state.temperature.fix(293.15)
model.flowsheet.pipe1.inlet_state.flow_mass.fix(550) #kg/s

#dummy objective
model.flowsheet.obj=pyo.Objective(expr=(model.flowsheet.pipe1.inlet_state.pressure-model.flowsheet.pipe1.outlet_state.pressure)**2)

#create compressor
model.flowsheet.compressor = idaesPressureChanger.PressureChanger(
        dynamic=False,
        property_package=model.flowsheet.paramBlock,
        compressor=True,
        thermodynamic_assumption=idaesPressureChanger.ThermodynamicAssumption.isentropic,
        )
model.flowsheet.compressor.display()
print(idaes.core.util.model_statistics.degrees_of_freedom(model))
#create constraints to connect pipe and compressor
model.flowsheet.comp_flow_con = pyo.Constraint(
        expr=model.flowsheet.compressor.inlet.flow_mass[0]==model.flowsheet.pipe1.outlet_state.flow_mass
        )
model.flowsheet.comp_inlet_pressure_con = pyo.Constraint(
        expr=model.flowsheet.compressor.inlet.pressure[0]==model.flowsheet.pipe1.outlet_state.pressure
        )
model.flowsheet.comp_temp_con = pyo.Constraint(
        expr=model.flowsheet.compressor.inlet.temperature[0]==model.flowsheet.pipe1.outlet_state.temperature
        )
#define other compressor constraints
model.flowsheet_comp_outlet_pressure_con = pyo.Constraint(
        expr=model.flowsheet.compressor.outlet.pressure[0]==350*100000 #35 MPa
        )
model.flowsheet.compressor.efficiency_isentropic.fix(0.85)

#add reservoir
model.flowsheet.reservoir = pyo.Block(rule=lambda b: reservoir_rule(b,reservoir_pressure=150*100000,reservoir_temperature=323,BOI=1.3,HCPV=1.5,PC_A=0.407215,PC_B=0.6972,GB_A=0.61805,GB_B=1.1294,GOR=0.51,paramBlock=model.flowsheet.paramBlock))
#connect to compressor
model.flowsheet.reservoir_inlet_flow_con = pyo.Constraint(
        expr=model.flowsheet.reservoir.inlet_state.flow_mass==model.flowsheet.compressor.outlet.flow_mass[0]
        )
model.flowsheet.reservoir_inlet_pressure_con = pyo.Constraint(
        expr=model.flowsheet.reservoir.inlet_state.pressure==model.flowsheet.compressor.outlet.pressure[0]
        )
model.flowsheet.reservoir_inlet_temp_con = pyo.Constraint(
        expr=model.flowsheet.reservoir.inlet_state.temperature==model.flowsheet.compressor.outlet.temperature[0]
        )


#manually scale
scale_mult = 0.1
#pipe scales
model.flowsheet.scaling_factor=pyo.Suffix(direction=pyo.Suffix.EXPORT)
model.flowsheet.scaling_factor[model.flowsheet.pipe1.diameter]=10*scale_mult
model.flowsheet.scaling_factor[model.flowsheet.pipe1.length]=1e-4*scale_mult
model.flowsheet.scaling_factor[model.flowsheet.pipe1.roughness]=1e2*scale_mult
model.flowsheet.scaling_factor[model.flowsheet.pipe1.area]=10*scale_mult
model.flowsheet.scaling_factor[model.flowsheet.pipe1.R]=1*scale_mult
model.flowsheet.scaling_factor[model.flowsheet.pipe1.inlet_velocity]=1*scale_mult
model.flowsheet.scaling_factor[model.flowsheet.pipe1.outlet_velocity]=1*scale_mult
model.flowsheet.scaling_factor[model.flowsheet.pipe1.average_velocity]=1*scale_mult
model.flowsheet.scaling_factor[model.flowsheet.pipe1.inlet_state.flow_mass]=1e-2*scale_mult
model.flowsheet.scaling_factor[model.flowsheet.pipe1.outlet_state.flow_mass]=1e-2*scale_mult
model.flowsheet.scaling_factor[model.flowsheet.pipe1.average_state.flow_mass]=1e-2*scale_mult
model.flowsheet.scaling_factor[model.flowsheet.pipe1.inlet_state.dens_mass]=1e-2*scale_mult
model.flowsheet.scaling_factor[model.flowsheet.pipe1.outlet_state.dens_mass]=1e-2*scale_mult
model.flowsheet.scaling_factor[model.flowsheet.pipe1.average_state.dens_mass]=1e-2*scale_mult
model.flowsheet.scaling_factor[model.flowsheet.pipe1.inlet_state.pressure]=1e-7*scale_mult
model.flowsheet.scaling_factor[model.flowsheet.pipe1.outlet_state.pressure]=1e-7*scale_mult
model.flowsheet.scaling_factor[model.flowsheet.pipe1.average_state.pressure]=1e-7*scale_mult
model.flowsheet.scaling_factor[model.flowsheet.pipe1.inlet_state.temperature]=1e-2*scale_mult
model.flowsheet.scaling_factor[model.flowsheet.pipe1.outlet_state.temperature]=1e-2*scale_mult
model.flowsheet.scaling_factor[model.flowsheet.pipe1.average_state.temperature]=1e-2*scale_mult
model.flowsheet.scaling_factor[model.flowsheet.pipe1.average_state.mw]=1e3*scale_mult
model.flowsheet.scaling_factor[model.flowsheet.pipe1.average_Z]=10*scale_mult
model.flowsheet.scaling_factor[model.flowsheet.pipe1.average_flux_mass]=1e-3*scale_mult
model.flowsheet.scaling_factor[model.flowsheet.pipe1.average_state.visc_d_phase["Vap"]]=1e6*scale_mult
model.flowsheet.scaling_factor[model.flowsheet.pipe1.average_reynolds]=1e-7*scale_mult
model.flowsheet.scaling_factor[model.flowsheet.pipe1.average_friction_factor]=100*scale_mult
model.flowsheet.scaling_factor[model.flowsheet.pipe1.isothermal]=1e-2*scale_mult
model.flowsheet.scaling_factor[model.flowsheet.pipe1.mass_bal]=1e4*scale_mult
model.flowsheet.scaling_factor[model.flowsheet.pipe1.hydraulic]=1e-3*scale_mult

#reservoir scales
#model.flowsheet.scaling_factor[model.flowsheet.reservoir.reservoir_pressure]=1e-7*scale_mult
#model.flowsheet.scaling_factor[model.flowsheet.reservoir.reservoir_temperature]=1e-2*scale_mult
model.flowsheet.scaling_factor[model.flowsheet.reservoir.BOI]=10*scale_mult
model.flowsheet.scaling_factor[model.flowsheet.reservoir.HCPV]=1e-2*scale_mult
model.flowsheet.scaling_factor[model.flowsheet.reservoir.PC_A]=0.1*scale_mult
model.flowsheet.scaling_factor[model.flowsheet.reservoir.PC_B]=0.1*scale_mult
model.flowsheet.scaling_factor[model.flowsheet.reservoir.GB_A]=1e3*scale_mult
model.flowsheet.scaling_factor[model.flowsheet.reservoir.GB_B]=1*scale_mult
model.flowsheet.scaling_factor[model.flowsheet.reservoir.GOR]=10*scale_mult
model.flowsheet.scaling_factor[model.flowsheet.reservoir.inlet_state.pressure]=1e-7*scale_mult
model.flowsheet.scaling_factor[model.flowsheet.reservoir.inlet_state.temperature]=1e-2*scale_mult
model.flowsheet.scaling_factor[model.flowsheet.reservoir.inlet_state.flow_mass]=1e-2*scale_mult
model.flowsheet.scaling_factor[model.flowsheet.reservoir.reservoir_state.pressure]=1e-7*scale_mult
model.flowsheet.scaling_factor[model.flowsheet.reservoir.reservoir_state.temperature]=1e-2*scale_mult
model.flowsheet.scaling_factor[model.flowsheet.reservoir.reservoir_state.flow_mass]=1e-2*scale_mult
model.flowsheet.scaling_factor[model.flowsheet.reservoir.mass_bal]=1e-2*scale_mult

model.flowsheet.scaling_factor[model.flowsheet.obj]=1e-13

#initialize compressor
model.flowsheet.compressor.initialize(outlvl=idaeslog.INFO)

scaled_model=pyo.TransformationFactory('core.scale_model').create_using(model)

solver = pyo.SolverFactory("ipopt")
solver.options['linear_solver']='ma27'
#solver.options['tol']=1e-6

result = solver.solve(scaled_model,logfile='ipopt_output.log',tee=True)

pyo.TransformationFactory('core.scale_model').propagate_solution(scaled_model,model)

result2=solver.solve(model,tee=True)

print(f'Pipe Results:')
print(f'pipe inlet velocity={pyo.value(model.flowsheet.pipe1.inlet_velocity)}')
print(f'pipe average mass density={pyo.value(model.flowsheet.pipe1.average_state.dens_mass)}')
print(f'pipe average Z={pyo.value(model.flowsheet.pipe1.average_Z)}')
print(f'mw={pyo.value(model.flowsheet.pipe1.average_state.mw)}')
print(f'pipe average velocity={pyo.value(model.flowsheet.pipe1.average_velocity)}')
print(f'pipe average reynolds={pyo.value(model.flowsheet.pipe1.average_reynolds)}')
print(f'pipe friction factor={pyo.value(model.flowsheet.pipe1.average_friction_factor)}')
print(f'pipe deltaP={pyo.value(model.flowsheet.pipe1.deltaP)}')
print()
print(f'Compressor Results')
model.flowsheet.compressor.report()
print()
print(f'reservoir results')
print(f'CO2 injection rate RB/s={pyo.value(model.flowsheet.reservoir.V_CO2_inj)}')
print(f'Oil production rate STB/s={pyo.value(model.flowsheet.reservoir.V_oil_prod)}')
print(f'Gas production rate RB/s={pyo.value(model.flowsheet.reservoir.V_gas_prod)}')
print(f'CO2 breakthrough rate RB/s={pyo.value(model.flowsheet.reservoir.V_CO2_brkth)}')
print(f'Total gas out RB/s={pyo.value(model.flowsheet.reservoir.V_gas_out)}')

#testing pipe model
"""
#only pipe, loop through inlet temperatures and get pressure drop
inlet_temps = [20,25,30,35,40,45,50,55,60]
model.flowsheet.pipe1.inlet_state.flow_mass.unfix()
model.flowsheet.velocity_constraint = pyo.Constraint(
        expr=model.flowsheet.pipe1.inlet_velocity==3
        )
deltaP_list = []
for inlet_temp in inlet_temps:
    model.flowsheet.pipe1.inlet_state.temperature.fix(inlet_temp+273.15)
    scaled_model = pyo.TransformationFactory('core.scale_model').create_using(model)
    result = solver.solve(scaled_model)
    pyo.TransformationFactory('core.scale_model').propagate_solution(scaled_model,model)
    result2=solver.solve(model)
    deltaP_list.append(pyo.value(model.flowsheet.pipe1.deltaP)/100000)
    print(f'inlet temp:{inlet_temp},deltaP:{deltaP_list[-1]}')
"""
