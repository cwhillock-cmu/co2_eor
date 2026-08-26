import pyomo.environ as pyo
import idaes.core
from pyomo.network import Arc
from pyomo.environ import units

from idaes.core import MomentumBalanceType
from co2_eor.splitter_unit import EnergySplittingType
from co2_eor.mixer_unit import MomentumMixingType

from co2_eor.SSLW import SSLWCosting, SSLWCostingData
from co2_eor.SSLW import CompressorType, CompressorDriveType, CompressorMaterial
from co2_eor.SSLW import VesselMaterial

from idaes.models.properties.modular_properties.base.generic_property import GenericParameterBlock

from co2_eor.MPF import configuration_vap_ideal, configuration_liq_ideal, configuration_VLE_ideal

from co2_eor.SSLW import SSLWCosting, SSLWCostingData

m = pyo.ConcreteModel()
m.fs = idaes.core.FlowsheetBlock(dynamic=False)

m.fs.ideal_vap = GenericParameterBlock(**configuration_vap_ideal)
m.fs.ideal_liq = GenericParameterBlock(**configuration_liq_ideal)
m.fs.ideal_VLe = GenericParameterBlock(**configuration_VLE_ideal)

m.fs.costing = SSLWCosting()

