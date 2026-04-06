import pyomo.environ as pyo
import idaes.core as idaescore
import pyomo.util as pyoutil
import idaes.models.properties.general_helmholtz as idaesHelmholtz
from co2_eor.mixer import MomentumMixingType
from co2_eor import mixer

#test block
m = pyo.ConcreteModel()
m.fs = idaescore.FlowsheetBlock(dynamic=False)
m.fs.props = idaesHelmholtz.HelmholtzParameterBlock(
        pure_component="CO2",amount_basis=idaesHelmholtz.AmountBasis.MASS,
        state_vars=idaesHelmholtz.StateVars.TPX,phase_presentation=idaesHelmholtz.PhaseType.L,
        #has_phase_equilibrium=False,
        )

m.fs.mixer = mixer(
        property_package=m.fs.props,
        inlet_list=["source1","source2"],
        momentum_mixing_type= MomentumMixingType.inequality
        )

m.fs.mixer.pprint()
