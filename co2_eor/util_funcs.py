#misc utility methods
import pyomo.environ as pyo
import pandas as pd

ipopt = pyo.SolverFactory('ipopt')
ipopt.options['linear_solver']='ma27'
ipopt.options['tol']=1E-8
ipopt.options['acceptable_tol']=1E-6
ipopt.options['halt_on_ampl_error']='yes'
bonmin = pyo.SolverFactory('bonmin')
bonmin.options['linear_solver']='ma97'
bonmin.options['tol']=1E-8
bonmin.options['acceptable_tol']=1E-6
bonmin.options['halt_on_ampl_error']='yes'
bonmin.options['bonmin.algorithm']='B-BB'
snopt = pyo.SolverFactory('snopt')
snopt.options['outlev']=2
snopt.options['major_iterations_limit'] = 10000
snopt.options['major_feasibility_tolerance']=1e-8
minos = pyo.SolverFactory('minos')
minos.options['outlev']=2
conopt=pyo.SolverFactory('conopt')
conopt.options['outlev']=3
conopt.options['logfreq']=1
conopt.options['hess']=0
conopt.options['iterlim']=500
baron = pyo.SolverFactory('baron')
LindoGlobal = pyo.SolverFactory('lindoglobal')
LGO = pyo.SolverFactory('LGO')
knitro = pyo.SolverFactory('knitro')
knitro.options['convex']=0
knitro.options['algorithm']=6
knitro.options['gradopt']=3
knitro.options['hessopt']=4
knitro.options['hessian_no_f']=0
knitro.options['honorbnds']=1
knitro.options['eval_fcga']=0
knitro.options['derivcheck']=0
knitro.options['cg_precond']=0
knitro.options['datacheck']=1

def export_compressor_df(compressor):
    """
    gets the important variables from an IDAES compressor object and returns a dataframe
    """
    data = {
        "flowrate (kg/s)":pyo.value(compressor.control_volume.properties_in[0].flow_mass),
        "inlet pressure (bar)":pyo.value(compressor.inlet.pressure[0]/100000),
        "inlet temperature (K)":pyo.value(compressor.control_volume.properties_in[0].temperature),
        "outlet pressure (bar)":pyo.value(compressor.outlet.pressure[0]/100000),
        "outlet temperature (K)":pyo.value(compressor.control_volume.properties_out[0].temperature),
        "Work Mechanical (kW)":pyo.value(compressor.work_mechanical[0]/1000),
    }
    return pd.DataFrame(data,index=[compressor.name])

def export_splitter_df(splitter):
    """
    gets the important variables from a splitter object and returns a dataframe
    """

    data = {}

    for j in splitter.mixed_state[0].component_list:
        data.update({f'inlet flowrate {j} (mol/s)':pyo.value(splitter.mixed_state[0].flow_mol_comp[j])})
    data.update({f'inlet pressure (bar)':pyo.value(splitter.mixed_state[0].pressure)/100000})
    data.update({f'inlet temperature (bar)':pyo.value(splitter.mixed_state[0].temperature)})
    for s in splitter.create_outlet_list():
        for j in getattr(splitter,s+"_state")[0].component_list:
            data.update({f'{s} flowrate {j} (mol/s)':pyo.value(getattr(splitter,s+"_state")[0].flow_mol_comp[j])})
        data.update({f'{s} pressure (bar)':pyo.value(getattr(splitter,s+"_state")[0].pressure)/100000})
        data.update({f'{s} temperature (bar)':pyo.value(getattr(splitter,s+"_state")[0].temperature)})

    return pd.DataFrame(data,index=[splitter.name])

def export_mixer_df(mixer):
    """
    gets the important variables from a mixer object and returns a dataframe
    """

    data = {}

    for s in mixer.create_inlet_list():
        for j in getattr(mixer,s+"_state")[0].component_list:
            data.update({f'{s} flowrate {j} (mol/s)':pyo.value(getattr(mixer,s+"_state")[0].flow_mol_comp[j])})
        data.update({f'{s} pressure (bar)':pyo.value(getattr(mixer,s+"_state")[0].pressure)/100000})
        data.update({f'{s} temperature (bar)':pyo.value(getattr(mixer,s+"_state")[0].temperature)})

    for j in mixer.mixed_state[0].component_list:
        data.update({f'inlet flowrate {j} (mol/s)':pyo.value(mixer.mixed_state[0].flow_mol_comp[j])})
    data.update({f'inlet pressure (bar)':pyo.value(mixer.mixed_state[0].pressure)/100000})
    data.update({f'inlet temperature (bar)':pyo.value(mixer.mixed_state[0].temperature)})

    return pd.DataFrame(data,index=[mixer.name])

def export_heater_df(heater):
    """
    gets the important variables from the IDAES heater unit with custom stuff added for co2_eor
    """
    data = {
        "flowrate (kg/s)":pyo.value(heater.control_volume.properties_in[0].flow_mass),
        "inlet pressure (bar)":pyo.value(heater.inlet.pressure[0]/100000),
        "inlet temperature (K)":pyo.value(heater.control_volume.properties_in[0].temperature),
        "outlet pressure (bar)":pyo.value(heater.outlet.pressure[0]/100000),
        "outlet temperature (K)":pyo.value(heater.control_volume.properties_out[0].temperature),
        "Heat Duty (kW)":pyo.value(heater.heat_duty[0]/1000),
        "Area (m2)":pyo.value(heater.area),
        "CW flowrate (m3/s)":pyo.value(heater.q_CW),
    }
    return pd.DataFrame(data,index=[heater.name])

        
def export_flowsheet_to_excel(flowsheet, filename):
    """
    Exports the data from pipeline, wellpad, and compressor units in the given
    IDAES flowsheetblock to an Excel file,
    along with overall flowsheet economics.
    """
    from idaes.models.unit_models.pressure_changer import PressureChanger
    from co2_eor.liq_pipe import liqPipe
    from co2_eor.gas_pipe import gasPipe 
    from co2_eor.processing_facility import processingFacility 
    from co2_eor.wellpattern_HH import wellpattern 
    from co2_eor.mpf_to_helmholtz import mpf_helmholtz_converter
    from co2_eor.mixer_unit import Mixer as mixer
    from co2_eor.splitter_unit import Separator as splitter
    from co2_eor.heater_unit import Heater as heater
    
    liqpipe_dfs = []
    gaspipe_dfs = []
    well_dfs = []
    comp_dfs = []
    processing_facility_dfs =[]
    mixer_splitter_dfs = []
    misc_dfs = []
    flowsheet_data = {}
    
    # Iterate over all unit model blocks in the flowsheet
    for blk_name, blk in flowsheet.component_map(pyo.Block).items():
        if isinstance(blk, liqPipe):
            liqpipe_dfs.append(blk.export_df())
            if hasattr(blk, "costing") and hasattr(blk.costing, "capital_cost"):
                flowsheet_data[f"{blk_name} capex"] = pyo.value(blk.costing.capital_cost)
        elif isinstance(blk, gasPipe):
            gaspipe_dfs.append(blk.export_df())
            if hasattr(blk, "costing") and hasattr(blk.costing, "capital_cost"):
                flowsheet_data[f"{blk_name} capex"] = pyo.value(blk.costing.capital_cost)
        elif isinstance(blk, wellpattern):
            well_dfs.append(blk.export_df())
            if hasattr(blk, "costing") and hasattr(blk.costing, "capital_cost"):
                flowsheet_data[f"{blk_name} capex"] = pyo.value(blk.costing.capital_cost)
        elif isinstance(blk, PressureChanger):
            comp_dfs.append(export_compressor_df(blk))
            if hasattr(blk, "costing") and hasattr(blk.costing, "capital_cost"):
                flowsheet_data[f"{blk_name} capex"] = pyo.value(blk.costing.capital_cost)
        elif isinstance(blk, processingFacility):
            processing_facility_dfs.append(blk.export_df())
            if hasattr(blk, "costing") and hasattr(blk.costing, "capital_cost"):
                flowsheet_data[f"{blk_name} capex"] = pyo.value(blk.costing.capital_cost)
        elif isinstance(blk, mpf_helmholtz_converter):
            misc_dfs.append(blk.export_df())
            if hasattr(blk, "costing") and hasattr(blk.costing, "capital_cost"):
                flowsheet_data[f"{blk_name} capex"] = pyo.value(blk.costing.capital_cost)
        elif isinstance(blk, heater):
            misc_dfs.append(export_heater_df(blk))
            if hasattr(blk, "costing") and hasattr(blk.costing, "capital_cost"):
                flowsheet_data[f"{blk_name} capex"] = pyo.value(blk.costing.capital_cost)
        elif isinstance(blk, mixer):
            mixer_splitter_dfs.append(export_mixer_df(blk))
        elif isinstance(blk, splitter):
            mixer_splitter_dfs.append(export_splitter_df(blk))
    
    # Economics
    if hasattr(flowsheet, "capex"):
        flowsheet_data["total capex"] = pyo.value(flowsheet.capex)
    if hasattr(flowsheet, "opex"):
        flowsheet_data["total opex"] = pyo.value(flowsheet.opex)
    if hasattr(flowsheet, "raw_mats"):
        flowsheet_data["total raw material cost"] = pyo.value(flowsheet.raw_mats)
    if hasattr(flowsheet, "revenue"):
        flowsheet_data["total revenue"] = pyo.value(flowsheet.revenue)
    if hasattr(flowsheet, "obj"):
        flowsheet_data["objective function"] = pyo.value(flowsheet.obj)
        
    flowsheet_df = pd.DataFrame(flowsheet_data, index=[1])
    
    with pd.ExcelWriter(filename) as writer:
        if liqpipe_dfs:
            pd.concat(liqpipe_dfs).to_excel(writer, sheet_name="LiquidPipeData", index=True)
        if gaspipe_dfs:
            pd.concat(gaspipe_dfs).to_excel(writer, sheet_name="GasPipeData", index=True)
        if well_dfs:
            pd.concat(well_dfs).to_excel(writer, sheet_name="WellData", index=True)
        if comp_dfs:
            pd.concat(comp_dfs).to_excel(writer, sheet_name="CompData", index=True)
        if processing_facility_dfs:
            pd.concat(processing_facility_dfs).to_excel(writer, sheet_name="ProcFacData", index=True)
        if misc_dfs:
            pd.concat(misc_dfs).to_excel(writer, sheet_name="miscUnitData", index=True)
        if mixer_splitter_dfs:
            pd.concat(mixer_splitter_dfs).to_excel(writer, sheet_name="MixSplits", index=True)
        flowsheet_df.to_excel(writer, sheet_name="FlowsheetData", index=True)

from idaes.core.util.initialization import propagate_state
from co2_eor.mixer_unit import Mixer as mixer
from co2_eor.splitter_unit import Separator as splitter
from co2_eor.wellpattern_HH import wellpattern
def custom_initialize_unit(unit):
    try:
        unit.initialize()
    except ValueError:
        #just in case, try deactivating unit feasibility problem
        try:
            unit.deactivate_feasibility_problem()
        except AttributeError:
            pass
        if isinstance(unit,(mixer, splitter, wellpattern)):
            #if unit is one of these, would take more effort to propagate the state, just return and hope for the best
            return
        print(f'call unit.initialize() failed, propagating state')
        #unit.display()
        #input()
        propagate_state(unit.inlet,unit.outlet)
 
#THIS FUNCTION IS AI-GENERATED 
from IPython.display import display 
def split_print_wide_df(df, max_cols=4, use_display=True):
    """
    Splits a wide pandas DataFrame into multiple narrower DataFrames 
    and automatically prints or displays each one.
    
    Args:
        df (pd.DataFrame): The input wide DataFrame.
        max_cols (int): The maximum number of columns per narrower DataFrame chunk.
        use_display (bool): If True, uses Jupyter's display() for rich HTML formatting. 
                            If False, uses standard text print().
    """
    # Extract the list of all columns
    cols = df.columns
    total_cols = len(cols)
    
    # Iterate through the columns in chunks of size `max_cols`
    for i in range(0, total_cols, max_cols):
        # Select the chunk of columns for this iteration
        chunk_cols = cols[i:i + max_cols]
        narrow_df = df[chunk_cols]
        
        # Print a header for context
        print(f"--- Columns {i+1} to {min(i+max_cols, total_cols)} ---")
        
        # Display or print the narrower DataFrame
        if use_display:
            try:
                # 'display' is the standard way to render DataFrames in Jupyter
                display(narrow_df)
            except NameError:
                # Fallback to standard print if display() isn't imported/available
                print(narrow_df)
        else:
            print(narrow_df)
            
        # Add a little spacing between the output dataframes
        print("\n")

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
def add_state_material_balances(self, balance_type=None, state_1=None, state_2=None,name='state_material_balances',doc='doc'):
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
    if hasattr(self, "state_material_balances"):
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

        setattr(self,name,state_material_balance_eq)
        


    elif balance_type == MaterialBalanceType.componentTotal:
        
        def state_material_balances_rule(t, j):
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