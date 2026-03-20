import pyomo.environ as pyo
from pyomo.common.config import ConfigBlock, ConfigValue, In, ListOf,Bool
from pyomo.network import Arc
from idaes.core import (
    ControlVolume0DBlock,
    declare_process_block_class,
    EnergyBalanceType,
    MomentumBalanceType,
    MaterialBalanceType,
    UnitModelBlockData,
    useDefault,
    FlowsheetBlock,
)
from idaes.core.util.config import is_physical_parameter_block
from idaes.core.util.scaling import set_scaling_factor
from idaes.core.scaling.autoscaling import AutoScaler
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.models.unit_models import Mixer, Separator, MomentumMixingType, MixingType, SplittingType,EnergySplittingType

#gather inputs global, for mixer, and for separator
def make_node_config_block(config):
    #global config values
    config.declare(
            "property_package",
            ConfigValue(default=useDefault,domain=is_physical_parameter_block)
            )
    config.declare(
            "property_package_args",
            ConfigBlock(implicit=True)
            )
    config.declare(
            "has_phase_equilibrium",
            ConfigValue(default=False,domain=Bool)
            )
    #config.declare(
    #        "dynamic",
    #        ConfigValue(domain=In([False]),default=True)
    #        )
    #config.declare(
    #        "has_holdup",
    #        ConfigValue(default=False,domain=In([False]))
    #        )
    #config values for mixer object
    config.declare(
            "inlet_list",
            ConfigValue(domain=ListOf(str))
            )
    config.declare(
            "num_inlets",
            ConfigValue(domain=int)
            )
    #config values for separator object
    config.declare(
            "outlet_list",
            ConfigValue(domain=ListOf(str))
            )
    config.declare(
            "num_outlets",
            ConfigValue(domain=int)
            )

#create a local flowsheet, add mixer and splitter units, and connect them to form node
def add_constraints(unit,name,config):
    unit.mixed_state = config.property_package.build_state_block(unit.flowsheet().time,defined_state=False)
    unit.mix = Mixer(
            property_package=config.property_package,
            inlet_list=config.inlet_list,
            momentum_mixing_type=MomentumMixingType.minimize,
            mixed_state_block=unit.mixed_state,
            )
    unit.split = Separator(
            property_package=config.property_package,
            ideal_separation=False,
            outlet_list=config.outlet_list,
            energy_split_basis=EnergySplittingType.equal_molar_enthalpy,
            momentum_balance_type=MomentumBalanceType.none,
            #momentum_balance_type=MomentumBalanceType.pressureTotal,
            mixed_state_block=unit.mixed_state,
            )
    #unit.arc = Arc(source=unit.mix.outlet,destination=unit.split.inlet)
    #pyo.TransformationFactory("network.expand_arcs").apply_to(unit)

    #splitter outlet pressure <= mixed state pressure
    @unit.split.Constraint(unit.flowsheet().time,unit.split.outlet_idx,)
    def pressure_inequality_constraint(b,t,o):
        o_block = getattr(unit.split,o+"_state")
        return unit.mixed_state[t].pressure>=o_block[t].pressure

#create node class
@declare_process_block_class("node")
class nodeData(UnitModelBlockData):
    CONFIG = UnitModelBlockData.CONFIG()
    make_node_config_block(CONFIG)

    def build(self):
        super(nodeData,self).build()
        add_constraints(self,"equations",self.config)

    def initialize(self):
        self.mix.initialize()
        self.split.initialize()






