#misc utility methods
import pyomo.environ as pyo
import pandas as pd

def export_compressor_df(compressor):
    """
    gets the important variables from an IDAES compressor object and returns a dataframe
    """
    data = {
        "flowrate (kg/s)":pyo.value(compressor.inlet.flow_mass[0]),
        "inlet pressure (bar)":pyo.value(compressor.inlet.pressure[0]/100000),
        "inlet temperature (K)":pyo.value(compressor.inlet.temperature[0]),
        "outlet pressure (bar)":pyo.value(compressor.outlet.pressure[0]/100000),
        "outlet temperature (K)":pyo.value(compressor.outlet.temperature[0]),
        "Work Mecahnical (kW)":pyo.value(compressor.work_mechanical[0]/1000),
    }
    return pd.DataFrame(data,index=[compressor.name])