#get data into dataframe
import pandas as pd
#df = pd.read_excel('temps/sensitivitytest3.xlsx',sheet_name=0)
df = pd.read_csv('temps/LHSdata1.csv')

clean_df = df[df['UOSSTAT2']==8]

print(clean_df)

input_labels = ['feed_co2_molflow',	'feed_ch4_molflow',	'feed_pressure',	'feed_temperature',	'recycle_solvent_molflow',
                	'HX1_outlet_temperature',	'flash3_outlet_pressure',	'column_stages',]
output_labels = ['recycle_outlet_co2_molflow',	'recycle_outlet_ch4_molflow',	'recycle_outlet_pressure',#	'recycle_outlet_temperature',
                 	'top_outlet_temperature',	'column_diameter',	'HX1_area',	'HX1_duty',	'pump1_volflow_in',	'pump1_work',	'comp1_work',
                        	'comp2_work',	'comp2_duty',	'flash1_volflow_in',	'flash2_volflow_in',	'flash3_volflow_in',	'makeup_solvent_molflow',]

surrogate_input_bounds = [
    [400,200,50e5,320,2000,260,6e5,10.5],
    [100,50,10e5,290,800,220,2e5,3.5]
]
input_bounds = {input_labels[i]: (surrogate_input_bounds[1][i],surrogate_input_bounds[0][i]) for i in range(len(input_labels))}

from idaes.core.surrogate.sampling.data_utils import split_training_validation
training_df, validation_df = split_training_validation(clean_df,0.8,seed=1)

import idaes.core
#"""
trainer = idaes.core.surrogate.AlamoTrainer(
    input_labels=input_labels,
    output_labels=output_labels,
    training_dataframe = training_df,
    input_bounds=input_bounds,
)
trainer.config.constant = 1
trainer.config.linfcns = 1
trainer.config.monomialpower = [2]
trainer.config.multi2power = [1,2]
trainer.config.logfcns = 1
trainer.config.expfcns = 1
trainer.config.ratiopower = [1]
trainer.config.screener = 0
trainer.config.modeler = 5
trainer.config.maxtime = 14400
trainer.config.maxterms = [67] * len(output_labels)
trainer.config.ZMIN = [0,0,0,0,0,0,0,0,0,0,0,-1e8,0,0,0,0]
trainer.config.xfactor = [1,1,1e5,1,1e3,1,1e5,1]
success, alm_surr, msg = trainer.train_surrogate()

surrogate_expressions = trainer._results['Model']
alm_surr = idaes.core.surrogate.AlamoSurrogate(surrogate_expressions, input_labels, output_labels,input_bounds)
model = alm_surr.save_to_file('temps/alamo_surrogate.json', overwrite=True)

from idaes.core.surrogate.plotting.sm_plotter import surrogate_scatter2D, surrogate_parity, surrogate_residual
#surrogate_scatter2D(alm_surr,validation_df,filename='temps/alamo_scatter.pdf',show=False)
surrogate_parity(alm_surr,validation_df,filename='temps/alamo_parity.pdf',show=False)
#surrogate_residual(alm_surr,validation_df,filename='temps/alamo_residual.pdf',show=False)
