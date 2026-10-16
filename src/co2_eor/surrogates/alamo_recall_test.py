import idaes.core
import pandas as pd

alm_surr = idaes.core.surrogate.AlamoSurrogate.load_from_file("temps/alamo_surrogate.json")

df = pd.read_excel("temps/LHS_data1_cut.xlsx")
input_df = df[alm_surr.input_labels()]
surrogate_output = alm_surr.evaluate_surrogate(input_df)

output_df = pd.concat([input_df,surrogate_output],axis=1)
output_df.to_excel("temps/LHSdata1_predicted.xlsx",index=False)