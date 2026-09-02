import pyomo.environ as pyo
import contextlib
from pyomo.environ import units
import idaes
import pyomo.util as pyoutil
import numpy as np
import idaes.models.properties.general_helmholtz as idaesHelmholtz
from idaes.models.properties.modular_properties.base.generic_property import GenericParameterBlock
import co2_eor.MPF.thermo_config as thermo_config
import idaes.logger as idaeslog
from co2_eor.util_funcs import ipopt

m = pyo.ConcreteModel()
m.fs = idaes.core.FlowsheetBlock(dynamic=False)
m.fs.props1 = GenericParameterBlock(**thermo_config.configuration_VLE_cubic)
m.fs.props2 = GenericParameterBlock(**thermo_config.configuration_vap_cubic)
m.fs.props3 = GenericParameterBlock(**thermo_config.configuration_liq_cubic)

m.fs.sb1 = m.fs.props3.build_state_block(has_phase_equilibrium=False,defined_state=True)

m.fs.sb1.flow_mol_comp['co2'].fix(100)
m.fs.sb1.flow_mol_comp['ch4'].fix(0)
m.fs.sb1.pressure.fix(550*100000)
m.fs.sb1.temperature.fix(300)
m.fs.sb1.test_expr = pyo.Expression(expr=m.fs.sb1.pressure-m.fs.sb1.pressure_crit)

with open('temps/test_liq_props_pprint.txt', 'w') as f:
    with contextlib.redirect_stdout(f):
        m.fs.sb1.pprint()

m.fs.sb1.initialize(outlvl=idaeslog.DEBUG)
res = ipopt.solve(m.fs.sb1,tee=True)
m.fs.sb1.display()
print(f'{pyo.value(m.fs.sb1.enth_mol)=}')
print(f'{pyo.value(m.fs.sb1.pressure_crit)=}')

input('paused')
m.fs.del_component(m.fs.sb1)

import pandas as pd

P_list = np.linspace(100,800,20)
f_ch4_list = np.flip(np.arange(0,105,5))

data_list=[]

for P in P_list:
    for f_ch4 in f_ch4_list:
        print(f'starting run P={P} bar, F_ch4 = {f_ch4} mol/s')
        m.fs.sb1 = m.fs.props1.build_state_block(has_phase_equilibrium=False,defined_state=True)
        m.fs.sb1.test_expr = pyo.Expression(expr=m.fs.sb1.pressure-m.fs.sb1.pressure_crit)
        m.fs.sb1.pressure.fix(P*100000)
        m.fs.sb1.flow_mol_comp['co2'].fix(100-f_ch4)
        m.fs.sb1.flow_mol_comp['ch4'].fix(f_ch4)
        try:
            m.fs.sb1.initialize(outlvl=idaeslog.WARNING)
            res = ipopt.solve(m.fs.sb1)
            is_optimal = (
                res.solver.status == pyo.SolverStatus.ok and 
                res.solver.termination_condition == pyo.TerminationCondition.optimal
            )
            assert is_optimal
            converged = True
            pressure_crit = pyo.value(m.fs.sb1.pressure_crit)
            enth_mol = pyo.value(m.fs.sb1.enth_mol)
        except (idaes.core.util.exceptions.InitializationError, AssertionError, ValueError) as e:
            converged = False
            pressure_crit = None
            enth_mol = None

        data_list.append({
            "P":P,
            'F_ch4':f_ch4,
            'F_co2':100-f_ch4,
            'converged':converged,
            'Critical P':pressure_crit,
            'Molar Enthalpy':enth_mol
            })
        m.fs.del_component(m.fs.sb1)

res_df = pd.DataFrame(data_list)
res_df.to_csv('temps/LiqPropsOutput.csv', index=False)