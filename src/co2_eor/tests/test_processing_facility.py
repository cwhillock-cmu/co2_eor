import pyomo.environ as pyo
import idaes
import pyomo.util as pyoutil
from pyomo.network import Arc
import idaes.models.properties.general_helmholtz as idaesHelmholtz
import idaes.models.unit_models.pressure_changer as idaesPressureChanger
import contextlib
from idaes.models.properties.modular_properties.base.generic_property import GenericParameterBlock
from co2_eor.MPF.thermo_config import configuration_vap
from co2_eor.processing_facility import processingFacility
from co2_eor.mpf_to_helmholtz import mpf_helmholtz_converter

m = pyo.ConcreteModel()
m.fs = idaes.core.FlowsheetBlock(dynamic=False)
m.fs.props1 = idaesHelmholtz.HelmholtzParameterBlock(
        pure_component="CO2",amount_basis=idaesHelmholtz.AmountBasis.MOLE,
        state_vars=idaesHelmholtz.StateVars.TPX,phase_presentation=idaesHelmholtz.PhaseType.MIX
        )
m.fs.props2 = GenericParameterBlock(**configuration_vap)

m.fs.processing_facility = processingFacility(
    property_package=m.fs.props2,
    has_phase_equilibrium=False,
    key_component='co2',
    recycle_pressure=75*100000,
)

m.fs.mpfConverter = mpf_helmholtz_converter(
    property_package_in=m.fs.props2,
    property_package_out=m.fs.props1,
    has_phase_equilibrium=False,
    conversion_type='convert_all_mol'
)
m.fs.s1 = Arc(source=m.fs.processing_facility.recycle,destination=m.fs.mpfConverter.inlet)

m.fs.comp1 = idaesPressureChanger.PressureChanger(
    property_package= m.fs.props1,
    dynamic=False,
    compressor=True,
    thermodynamic_assumption=idaesPressureChanger.ThermodynamicAssumption.isentropic,
)
m.fs.comp1.efficiency_isentropic.fix(0.85)

m.fs.s2 = Arc(source=m.fs.mpfConverter.outlet,destination=m.fs.comp1.inlet)

pyo.TransformationFactory("network.expand_arcs").apply_to(m)

#m.fs.processing_facility.inlet.flow_mol[0].fix(9)
#m.fs.processing_facility.inlet.mol_frac_comp[0,'co2'].fix(0.9)
#m.fs.processing_facility.inlet.mol_frac_comp[0,'ch4'].fix(0.1)
m.fs.processing_facility.inlet.flow_mol_comp[0,'co2'].fix(2.0967739836727803)
m.fs.processing_facility.inlet.flow_mol_comp[0,'ch4'].fix(0.0267751879149411)
m.fs.processing_facility.inlet.pressure[0].fix(9*100000)
m.fs.processing_facility.inlet.temperature[0].fix(302)

m.fs.comp1.outlet.pressure[0].fix(100*100000)

with open('temps/test_processing_facility_pprint.txt', 'w') as f:
    with contextlib.redirect_stdout(f):
        m.pprint()

print(f'D.o.F ={idaes.core.util.model_statistics.degrees_of_freedom(m)}')
print(f'processing facility D.o.F ={idaes.core.util.model_statistics.degrees_of_freedom(m.fs.processing_facility)}')
print(f'converter D.o.F ={idaes.core.util.model_statistics.degrees_of_freedom(m.fs.mpfConverter)}')

m.fs.processing_facility.control_volume.properties_in.initialize()
m.fs.processing_facility.initialize(display_after=False)
m.fs.mpfConverter.initialize()

from co2_eor.util_funcs import ipopt
ripopt = pyo.SolverFactory('ripopt')
#scale model
scaled_m = pyo.TransformationFactory("core.scale_model").create_using(m)
#solve flowsheet
res=ipopt.solve(scaled_m,tee=True,logfile='ipopt_output.log')
#unscale model
pyo.TransformationFactory("core.scale_model").propagate_solution(scaled_m,m)

with open('temps/test_processing_facility_postsolve_display.txt', 'w') as f:
    with contextlib.redirect_stdout(f):
        m.display()

from co2_eor.util_funcs import split_print_wide_df, export_compressor_df
split_print_wide_df(m.fs.processing_facility.export_df(),max_cols=4)
split_print_wide_df(export_compressor_df(m.fs.comp1),max_cols=4)