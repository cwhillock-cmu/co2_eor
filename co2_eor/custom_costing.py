import pyomo.environ as pyo
from idaes.core import ProcessBlockData, declare_process_block_class
from pyomo.common.config import ConfigBlock, ConfigValue

def make_config_block(config):
    config.declare(
        "costing_method",
        ConfigValue(),
    )
    config.declare(
        "costing_method_args",
        ConfigValue(),
    )

@declare_process_block_class("customCostingBlock")
class customCostingBlockData(ProcessBlockData):
    config = ProcessBlockData.CONFIG()
    make_config_block(config)

    def build(self):
        super(customCostingBlockData,self).build()
        method = self.config.costing_method

        method(self,**self.config.costing_method_args)

def cost_processing_facility_simple(blk,F_in,F_recycle,oil_price):
    #simplified processing facility economics
    capex_coefficient = pyo.Param(initialize=42.35) #$/(tonne/yr)

    blk.capital_cost = pyo.Expression(
        expr=capex_coefficient*(F_recycle*3600*24*365/1000)
    )

    blk.opex = pyo.Expression(
        expr=oil_price*0.01*F_in/1.98*35.3147/1000 #kg/s -> Sm3/s -> Sft3/s -> MSCF/s
    )



    
