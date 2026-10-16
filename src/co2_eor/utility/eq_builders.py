import pyomo.environ as pyo

def add_state_material_balances(self, balance_type=None, state_1=None, state_2=None,name='state_material_balances',doc='doc'):
    from idaes.core.base.property_base import StateBlock
    from idaes.core.util.exceptions import (
        BurntToast,
        ConfigurationError,
        BalanceTypeNotSupportedError,
        InitializationError,
    )
    from idaes.core.base.control_volume_base import (
        ControlVolumeBlockData,
        FlowDirection,
        MaterialBalanceType,
    ) 
    """
    Method to add material balances linking two State Blocks in a Unit
    Model. This method is not intended to replace Control Volumes, but
    to automate writing material balances linking isolated State Blocks
    in those models where this is required.

    Args:
        balance_type - a MaterialBalanceType Enum indicating the type
                        of material balances to write
        state_1 - first State Block to be linked by balances
        state_2 - second State Block to be linked by balances

    Returns:
        None

    this method is modified to create the constraint and add it to the unit instead of using decorators.
    allows for a custom name.
    uses the self keyword for the unit because of convenience - not a class function.
    """
    # Confirm that both state blocks stem from the same parameter block
    if not isinstance(state_1, StateBlock):
        raise ConfigurationError(
            "{} state_1 argument to add_state_material_balances "
            "was not an instance of a State Block.".format(self.name)
        )

    if not isinstance(state_2, StateBlock):
        raise ConfigurationError(
            "{} state_2 argument to add_state_material_balances "
            "was not an instance of a State Block.".format(self.name)
        )

    # Check that no constraint with the same name exists
    # We will only support using this method once per Block
    if hasattr(self, name):
        raise AttributeError(
            "{} a set of constraints named state_material_balances "
            "already exists in the current UnitModel. To avoid "
            "confusion, add_state_material_balances is only supported "
            "once per UnitModel.".format(self.name)
        )

    # Get a representative time point for testing
    rep_time = self.flowsheet().time.first()
    if state_1[rep_time].params is not state_2[rep_time].params:
        raise ConfigurationError(
            "{} add_state_material_balances method was provided with "
            "State Blocks are not linked to the same "
            "instance of a Physical Parameter Block. This method "
            "only supports linking State Blocks from the same "
            "Physical Parameter Block.".format(self.name)
        )

    if balance_type == MaterialBalanceType.useDefault:
        balance_type = state_1[rep_time].default_material_balance_type()

    phase_list = state_1.phase_list
    component_list = state_1.component_list
    pc_set = state_1.phase_component_set

    if balance_type == MaterialBalanceType.componentPhase:
        # TODO : Should we include an optional phase equilibrium term here
        # to allow for systems where a phase-transition may occur?
        
        def state_material_balances_rule(b, t, p, j):
            return state_1[t].get_material_flow_terms(p, j) == state_2[
                t
            ].get_material_flow_terms(p, j)
        
        state_material_balance_eq = pyo.Constraint(
            self.flowsheet().time,
            pc_set,
            doc=doc,
            rule=state_material_balances_rule
        )

    elif balance_type == MaterialBalanceType.componentTotal:
        
        def state_material_balances_rule(b, t, j):
            return sum(
                state_1[t].get_material_flow_terms(p, j)
                for p in phase_list
                if (p, j) in pc_set
            ) == sum(
                state_2[t].get_material_flow_terms(p, j)
                for p in phase_list
                if (p, j) in pc_set
            )
    
        state_material_balance_eq = pyo.Constraint(
            self.flowsheet().time,
            component_list,
            doc=doc,
            rule=state_material_balances_rule
        )
        

    elif balance_type == MaterialBalanceType.total:

        def state_material_balances_rule(t):
            return sum(
                state_1[t].get_material_flow_terms(p, j) for p, j in pc_set
            ) == sum(state_2[t].get_material_flow_terms(p, j) for p, j in pc_set)
        
        state_material_balance_eq = pyo.Constraint(
            self.flowsheet().time,
            doc=doc,
            rule=state_material_balances_rule
        )

    elif balance_type == MaterialBalanceType.elementTotal:
        raise BalanceTypeNotSupportedError(
            "{} add_state_material_balances does not support "
            "MaterialBalanceType.elementTotal.".format(self.name)
        )
    elif balance_type == MaterialBalanceType.none:
        raise BalanceTypeNotSupportedError(
            "{} add_state_material_balances does not support "
            "MaterialBalanceType.None.".format(self.name)
        )
    else:
        raise BurntToast(
            "{} add_state_material_balances received an unexpected "
            "argument for balance_type. This should never happen. Please "
            "contact the IDAES developers with this bug.".format(self.name)
        )
    
    setattr(self,name,state_material_balance_eq)

def add_pressure_thickness_expression(blk,pressure=None,diameter=None,allowable_stress=None,overdesign=1.1):
    pressure_thickness_expression = pyo.Expression(expr=overdesign*pressure*diameter/(2*(allowable_stress-overdesign*pressure)))
    setattr(blk,"thickness",pressure_thickness_expression)