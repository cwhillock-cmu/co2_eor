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

m = pyo.ConcreteModel()
m.fs = idaes.core.FlowsheetBlock(dynamic=False)
m.fs.props1 = GenericParameterBlock(**thermo_config.configuration_VLE_cubic)
m.fs.props2 = GenericParameterBlock(**thermo_config.configuration_vap_cubic)
m.fs.props3 = GenericParameterBlock(**thermo_config.configuration_liq_cubic)

m.fs.sb1 = m.fs.props1.build_state_block(has_phase_equilibrium=True,defined_state=True)

m.fs.sb1.flow_mol_comp['co2'].fix(100)
m.fs.sb1.flow_mol_comp['ch4'].fix(20)
m.fs.sb1.pressure.fix(50*100000)
m.fs.sb1.temperature.fix(300)

with open('temps/test_VLE_pprint.txt', 'w') as f:
    with contextlib.redirect_stdout(f):
        m.fs.sb1.pprint()

m.fs.sb1.initialize(outlvl=idaeslog.DEBUG)

m.fs.sb1.display()

input('paused')
m.fs.del_component(m.fs.sb1)

import pandas as pd
from co2_eor.util_funcs import ipopt
P_list = np.linspace(20,160,15)
f_ch4_list = np.flip(np.arange(5,100,5))

data_list=[]

for P in P_list:
    for f_ch4 in f_ch4_list:
        print(f'starting run P={P} bar, F_ch4 = {f_ch4} mol/s')
        m.fs.sb1 = m.fs.props1.build_state_block(has_phase_equilibrium=True,defined_state=True)
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
            phase_frac_vap = pyo.value(m.fs.sb1.phase_frac['Vap'])
            phase_frac_liq = pyo.value(m.fs.sb1.phase_frac['Liq'])
            T_bub = pyo.value(m.fs.sb1.temperature_bubble[('Vap','Liq')])
            T_dew = pyo.value(m.fs.sb1.temperature_dew[('Vap','Liq')])
            T_eq = pyo.value(m.fs.sb1._teq[('Vap','Liq')])

        except (idaes.core.util.exceptions.InitializationError, AssertionError, ValueError) as e:
            phase_frac_vap = None
            phase_frac_liq = None
            T_bub = None
            T_dew = None
            T_eq = None

        data_list.append({
            "P":P,
            'F_ch4':f_ch4,
            'F_co2':100-f_ch4,
            'Phase Frac Vap':phase_frac_vap,
            'Phase Frac Liq':phase_frac_liq,
            'T Bubble':T_bub,
            'T dew':T_dew,
            'T eq':T_eq
            })
        m.fs.del_component(m.fs.sb1)

res_df = pd.DataFrame(data_list)
res_df.to_csv('temps/LogBubbleDewoutput.csv', index=False)