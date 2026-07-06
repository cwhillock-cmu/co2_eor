import pyomo.environ as pyo
import idaes
import pyomo.util as pyoutil
from pyomo.network import Arc
import idaes.models.properties.general_helmholtz as idaesHelmholtz
import contextlib
from idaes.models.properties.modular_properties.base.generic_property import GenericParameterBlock
from co2_eor.wellpattern_r3 import wellpattern
from co2_eor.MPF.thermo_config import configuration
from co2_eor.gas_pipe import gasPipe

m = pyo.ConcreteModel()
m.fs = idaes.core.FlowsheetBlock(dynamic=False)
m.fs.props1 = idaesHelmholtz.HelmholtzParameterBlock(
        pure_component="CO2",amount_basis=idaesHelmholtz.AmountBasis.MASS,
        state_vars=idaesHelmholtz.StateVars.TPX,phase_presentation=idaesHelmholtz.PhaseType.L
        )
m.fs.props2 = GenericParameterBlock(**configuration)

m.fs.wellpattern = wellpattern(
        primary_property_package=m.fs.props1,
        secondary_property_package=m.fs.props2,
        temperature=305,
        pressure=138*100000,
        Boi=1.22,
        GOR=4.42,
        PC_A=1.386,
        PC_B=2.507,
        GB_A=0.710376,
        GB_B=1.8315,
        kovr=6.9e-13,
        use_correction_factor=False,
        SC_B=789,
        IR_base=0.001968,
        depth=1000,
        outlet_pressure=30*100000,
        outlet_temperature=300,
        )

m.fs.pipe = gasPipe(
        property_package=m.fs.props2,
        length=5000,
        average_pressure_type='nonlinear',
)
m.fs.pipe.diameter.fix(0.5)
m.fs.pipe.roughness.fix(0.0475e-3)

m.fs.s1 = Arc(source=m.fs.wellpattern.outlet,destination=m.fs.pipe.inlet)

pyo.TransformationFactory("network.expand_arcs").apply_to(m)

with open('temps/test_wellpatternv3_gas_pipe_pprint.txt', 'w') as f:
    with contextlib.redirect_stdout(f):
        m.pprint()

m.fs.wellpattern.inlet.temperature.fix(300)
m.fs.wellpattern.inlet.pressure.fix(100*100000)
m.fs.wellpattern.HCPV.fix(0.5)

print(f'D.o.F ={idaes.core.util.model_statistics.degrees_of_freedom(m)}')

m.fs.wellpattern.initialize()
m.fs.pipe.initialize()
from co2_eor.util_funcs import ipopt
#scale model
scaled_m = pyo.TransformationFactory("core.scale_model").create_using(m)
#solve flowsheet
res=ipopt.solve(scaled_m,tee=True,logfile='ipopt_output.log')
#unscale model
pyo.TransformationFactory("core.scale_model").propagate_solution(scaled_m,m)

with open('temps/test_wellpatternv3_gas_pipe_display.txt', 'w') as f:
    with contextlib.redirect_stdout(f):
        m.display()

from co2_eor.util_funcs import split_print_wide_df
split_print_wide_df(m.fs.wellpattern.export_df())
split_print_wide_df(m.fs.pipe.export_df())