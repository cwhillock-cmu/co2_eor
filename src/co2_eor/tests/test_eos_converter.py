import pyomo.environ as pyo
import idaes
import pyomo.util as pyoutil
from pyomo.network import Arc
import idaes.models.properties.general_helmholtz as idaesHelmholtz
import contextlib
from idaes.models.properties.modular_properties.base.generic_property import GenericParameterBlock
from co2_eor.MPF.thermo_config import configuration
from co2_eor.processing_facility import processingFacility
from co2_eor.mpf_to_helmholtz import mpf_helmholtz_converter

m = pyo.ConcreteModel()
m.fs = idaes.core.FlowsheetBlock(dynamic=False)
m.fs.props1 = idaesHelmholtz.HelmholtzParameterBlock(
        pure_component="CO2",amount_basis=idaesHelmholtz.AmountBasis.MASS,
        state_vars=idaesHelmholtz.StateVars.TPX,phase_presentation=idaesHelmholtz.PhaseType.L
        )
m.fs.props2 = GenericParameterBlock(**configuration)

m.fs.mpfConverter = mpf_helmholtz_converter(
    property_package_in=m.fs.props2,
    property_package_out=m.fs.props1,
    has_phase_equilibrium=False,
    conversion_type='convert_co2'
)

m.fs.mpfConverter.inlet.flow_mass[0].fix(9)
m.fs.mpfConverter.inlet.mass_frac_comp[0,'co2'].fix(0.9)
m.fs.mpfConverter.inlet.pressure[0].fix(2*100000)
m.fs.mpfConverter.inlet.temperature[0].fix(300)

with open('temps/test_eos_conversion_pprint.txt', 'w') as f:
    with contextlib.redirect_stdout(f):
        m.pprint()

print(f'D.o.F ={idaes.core.util.model_statistics.degrees_of_freedom(m)}')

m.fs.mpfConverter.initialize(display_after=True)