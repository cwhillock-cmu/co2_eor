import pyomo.environ as pyo
import idaes
import pyomo.util as pyoutil
from pyomo.network import Arc
import idaes.models.properties.general_helmholtz as idaesHelmholtz
import contextlib
from idaes.models.properties.modular_properties.base.generic_property import GenericParameterBlock
from co2_eor.thermo_config import configuration
from co2_eor.processing_facility import processingFacility
from co2_eor.mpf_to_helmholtz import mpf_helmholtz_converter

m = pyo.ConcreteModel()
m.fs = idaes.core.FlowsheetBlock(dynamic=False)
m.fs.props1 = idaesHelmholtz.HelmholtzParameterBlock(
        pure_component="CO2",amount_basis=idaesHelmholtz.AmountBasis.MASS,
        state_vars=idaesHelmholtz.StateVars.TPX,phase_presentation=idaesHelmholtz.PhaseType.L
        )
m.fs.props2 = GenericParameterBlock(**configuration)

m.fs.processing_facility = processingFacility(
    property_package=m.fs.props2,
    has_phase_equilibrium=False,
    key_component='co2',
)

m.fs.mpfConverter = mpf_helmholtz_converter(
    property_package_in=m.fs.props2,
    property_package_out=m.fs.props1,
    has_phase_equilibrium=False,
    conversion_type='convert_co2'
)
m.fs.s1 = Arc(source=m.fs.processing_facility.recycle,destination=m.fs.mpfConverter.inlet)
pyo.TransformationFactory("network.expand_arcs").apply_to(m)

m.fs.processing_facility.inlet.flow_mass[0].fix(9)
m.fs.processing_facility.inlet.mass_frac_comp[0,'co2'].fix(0.9)
m.fs.processing_facility.inlet.mass_frac_comp[0,'ch4'].fix(0.1)
m.fs.processing_facility.inlet.pressure[0].fix(9*100000)
m.fs.processing_facility.inlet.temperature[0].fix(302)

with open('temps/test_processing_facility_pprint.txt', 'w') as f:
    with contextlib.redirect_stdout(f):
        m.pprint()

print(f'D.o.F ={idaes.core.util.model_statistics.degrees_of_freedom(m)}')
print(f'processing facility D.o.F ={idaes.core.util.model_statistics.degrees_of_freedom(m.fs.processing_facility)}')
print(f'converter D.o.F ={idaes.core.util.model_statistics.degrees_of_freedom(m.fs.mpfConverter)}')

m.fs.processing_facility.initialize(display_after=True)
m.fs.mpfConverter.initialize()

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
split_print_wide_df(m.fs.processing_facility.export_df())