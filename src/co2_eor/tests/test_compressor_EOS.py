import pyomo.environ as pyo
from pyomo.environ import units
import idaes.core as idaescore
import pyomo.util as pyoutil
import numpy as np
import idaes.models.properties.general_helmholtz as idaesHelmholtz
from idaes.models.properties.modular_properties.base.generic_property import GenericParameterBlock
from co2_eor.MPF.thermo_config import configuration_VLE

m = pyo.ConcreteModel()
m.fs = idaescore.FlowsheetBlock(dynamic=False)
m.fs.props1 = idaesHelmholtz.HelmholtzParameterBlock(
        pure_component="CO2",amount_basis=idaesHelmholtz.AmountBasis.MOLE,
        state_vars=idaesHelmholtz.StateVars.PH,phase_presentation=idaesHelmholtz.PhaseType.MIX,
        #has_phase_equilibrium=False,
        )

import idaes.models.unit_models.pressure_changer as idaesPressureChanger

m.fs.comp1 = idaesPressureChanger.PressureChanger(
    property_package= m.fs.props1,
    dynamic=False,
    compressor=True,
    thermodynamic_assumption=idaesPressureChanger.ThermodynamicAssumption.isentropic,
)
m.fs.comp1.efficiency_isentropic.fix(0.85)

m.fs.comp1.inlet.pressure[0].fix(5*100000)
#m.fs.comp1.inlet.temperature[0].fix(298)
m.fs.comp1.inlet.enth_mol[0].fix(m.fs.props1.htpx(T=298*units.K,p=5*100000*units.Pa,amount_basis=idaesHelmholtz.AmountBasis.MOLE))
m.fs.comp1.inlet.flow_mol[0].fix(500)
m.fs.comp1.outlet.pressure[0].fix(100*100000)
#m.fs.comp1.inlet.vapor_frac[0].fix(1)

print(f'DoF={idaescore.util.model_statistics.degrees_of_freedom(m)}')

solver = pyo.SolverFactory('ipopt')

m.fs.comp1.initialize()

#scale model
scaled_m = pyo.TransformationFactory("core.scale_model").create_using(m)
#solve flowsheet
res=solver.solve(scaled_m,tee=True,logfile='ipopt_output.log')
#unscale model
pyo.TransformationFactory("core.scale_model").propagate_solution(scaled_m,m)

m.display()
print(f'inlet T: {pyo.value(m.fs.comp1.control_volume.properties_in[0].temperature)}')
print(f'outlet T: {pyo.value(m.fs.comp1.control_volume.properties_out[0].temperature)}')
print(f'inlet V frac: {pyo.value(m.fs.comp1.control_volume.properties_in[0].vapor_frac)}')
print(f'outlet V frac: {pyo.value(m.fs.comp1.control_volume.properties_out[0].vapor_frac)}')