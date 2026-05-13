#misc utility methods
import pyomo.environ as pyo
import pandas as pd

ipopt = pyo.SolverFactory('ipopt')
ipopt.options['linear_solver']='ma97'
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
conopt.options['outlev']=2
conopt.options['logfreq']=1
conopt.options['hess']=0
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
        "flowrate (kg/s)":pyo.value(compressor.inlet.flow_mass[0]),
        "inlet pressure (bar)":pyo.value(compressor.inlet.pressure[0]/100000),
        "inlet temperature (K)":pyo.value(compressor.control_volume.properties_in[0].temperature),
        "outlet pressure (bar)":pyo.value(compressor.outlet.pressure[0]/100000),
        "outlet temperature (K)":pyo.value(compressor.control_volume.properties_out[0].temperature),
        "Work Mecahnical (kW)":pyo.value(compressor.work_mechanical[0]/1000),
    }
    return pd.DataFrame(data,index=[compressor.name])

def export_flowsheet_to_excel(flowsheet, filename):
    """
    Exports the data from pipeline, wellpad, and compressor units in the given
    IDAES flowsheetblock to an Excel file,
    along with overall flowsheet economics.
    """
    from co2_eor.pipeline import pipeline
    from co2_eor.wellpad import wellpad
    from idaes.models.unit_models.pressure_changer import PressureChanger
    
    pipe_dfs = []
    well_dfs = []
    comp_dfs = []
    flowsheet_data = {}
    
    # Iterate over all unit model blocks in the flowsheet
    for blk_name, blk in flowsheet.component_map(pyo.Block).items():
        if isinstance(blk, pipeline):
            pipe_dfs.append(blk.export_df())
            if hasattr(blk, "costing") and hasattr(blk.costing, "capital_cost"):
                flowsheet_data[f"{blk_name} capex"] = pyo.value(blk.costing.capital_cost)
        elif isinstance(blk, wellpad):
            well_dfs.append(blk.export_df())
            if hasattr(blk, "costing") and hasattr(blk.costing, "capital_cost"):
                flowsheet_data[f"{blk_name} capex"] = pyo.value(blk.costing.capital_cost)
        elif isinstance(blk, PressureChanger):
            comp_dfs.append(export_compressor_df(blk))
            if hasattr(blk, "costing") and hasattr(blk.costing, "capital_cost"):
                flowsheet_data[f"{blk_name} capex"] = pyo.value(blk.costing.capital_cost)
    
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
        if pipe_dfs:
            pd.concat(pipe_dfs).to_excel(writer, sheet_name="PipeData", index=True)
        if well_dfs:
            pd.concat(well_dfs).to_excel(writer, sheet_name="WellData", index=True)
        if comp_dfs:
            pd.concat(comp_dfs).to_excel(writer, sheet_name="CompData", index=True)
        flowsheet_df.to_excel(writer, sheet_name="FlowsheetData", index=True)

from idaes.core.util.initialization import propagate_state
from co2_eor.pipeline import pipeline
from co2_eor.wellpad import wellpad
def fs_initializer_function(unit):
    print(f'initializing {unit.name}')
    if isinstance(unit, pipeline):
        try:
            res = unit.initialize()
            pyo.assert_optimal_termination(res)
        except (ValueError, RuntimeError):
            unit.deactivate_feasibility_problem()
            print(f'solve failed, propagating state')
            propagate_state(unit.inlet,unit.outlet)
    elif isinstance(unit,wellpad):
        try:
            res = unit.initialize()
            pyo.assert_optimal_termination(res)
        except (ValueError, RuntimeError):
            unit.deactivate_feasibility_problem()
            print(f'solve failed, propagating state')
            propagate_state(unit.inlet,unit.outlet)
    else:
        try:
            unit.initialize()
        except ValueError:
            print(f'solve failed, propagating state')
            propagate_state(unit.inlet,unit.outlet)