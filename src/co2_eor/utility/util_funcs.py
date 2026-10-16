#misc utility methods
import pyomo.environ as pyo
import pandas as pd

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
    from co2_eor.wellpattern_unit import wellpattern 
    from co2_eor.mpf_to_helmholtz import mpf_helmholtz_converter
    from co2_eor.mixer_unit import Mixer as mixer
    from co2_eor.splitter_unit import Separator as splitter
    from co2_eor.heater_unit import Heater as heater
    from co2_eor.surrogates.rectisol_surrogate import rectisolSurrogate
    
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
        elif isinstance(blk, rectisolSurrogate):
            processing_facility_dfs.append(blk.export_df())
            if hasattr(blk, "costing") and hasattr(blk.costing, "total_bare_module_cost"):
                flowsheet_data[f"{blk_name} TBMC"] = pyo.value(blk.costing.total_bare_module_cost)
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
