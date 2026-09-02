import pyomo.environ as pyo
import contextlib
from pyomo.environ import units
import idaes.core as idaescore
import pyomo.util as pyoutil
import numpy as np
import idaes.models.properties.general_helmholtz as idaesHelmholtz
from idaes.models.properties.modular_properties.base.generic_property import GenericParameterBlock
import co2_eor.MPF.thermo_config as thermo_config

m = pyo.ConcreteModel()
m.fs = idaescore.FlowsheetBlock(dynamic=False)
m.fs.props1 = GenericParameterBlock(**thermo_config.configuration_VLE_cubic)
m.fs.props2 = GenericParameterBlock(**thermo_config.configuration_vap_cubic)
m.fs.props3 = GenericParameterBlock(**thermo_config.configuration_liq_cubic)

import idaes.models.unit_models.pressure_changer as idaesPressureChanger

m.fs.comp1 = idaesPressureChanger.PressureChanger(
    property_package= m.fs.props1,
    dynamic=False,
    compressor=True,
    thermodynamic_assumption=idaesPressureChanger.ThermodynamicAssumption.isentropic,
    has_phase_equilibrium=False,
)
m.fs.comp1.efficiency_isentropic.fix(0.85)

m.fs.comp1.inlet.pressure[0].fix(200*100000)
m.fs.comp1.inlet.temperature[0].fix(300)
#m.fs.comp1.inlet.enth_mol[0].fix(m.fs.props1.htpx(T=298*units.K,p=5*100000*units.Pa,amount_basis=idaesHelmholtz.AmountBasis.MOLE))
m.fs.comp1.inlet.flow_mol_comp[0,'co2'].fix(100)
m.fs.comp1.inlet.flow_mol_comp[0,'ch4'].fix(0)
m.fs.comp1.outlet.pressure[0].fix(550*100000)

print(f'DoF={idaescore.util.model_statistics.degrees_of_freedom(m)}')

from co2_eor.util_funcs import ipopt, conopt
import idaes.logger as idaeslog
m.fs.comp1.initialize(outlvl=idaeslog.DEBUG)

#scale model
scaled_m = pyo.TransformationFactory("core.scale_model").create_using(m)
#solve flowsheet
res=ipopt.solve(scaled_m,tee=True)
#unscale model
pyo.TransformationFactory("core.scale_model").propagate_solution(scaled_m,m)

with open('temps/test_compressor_VLE_display.txt', 'w') as f:
    with contextlib.redirect_stdout(f):
        m.display()