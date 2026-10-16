import idaes.core
import pandas as pd
df = pd.read_excel('temps/LHSdata123.xlsx')

clean_df = df[df['UOSSTAT2']==8]

print(clean_df)

input_labels = ['feed_co2_molflow',	'feed_ch4_molflow',	'feed_pressure',	'feed_temperature',	'recycle_solvent_molflow',
                	'HX1_outlet_temperature',	'flash3_outlet_pressure',	'column_stages',]
output_labels = ['recycle_outlet_co2_molflow',	'recycle_outlet_ch4_molflow',	'recycle_outlet_pressure',	'recycle_outlet_temperature',
                 	'top_outlet_temperature',	'column_diameter',	'HX1_area',	'HX1_duty',	'pump1_volflow_in',	'pump1_work',	'comp1_work',
                        	'comp2_work',	'comp2_duty',	'flash1_volflow_in',	'flash2_volflow_in',	'flash3_volflow_in',	'makeup_solvent_molflow',]

surrogate_input_bounds = [
    [400,200,50e5,320,2500,260,10e5,10.5],
    [100,50,10e5,290,800,220,2e5,3.5]
]

input_bounds = {input_labels[i]: (surrogate_input_bounds[1][i],surrogate_input_bounds[0][i]) for i in range(len(input_labels))}

from idaes.core.surrogate.sampling.data_utils import split_training_validation
training_df, validation_df = split_training_validation(clean_df,0.70,seed=1)

input_data = training_df.iloc[:, 2:10]
output_data = training_df.iloc[:, 10:27]

import tensorflow as tf

# selected settings for regression
activation, optimizer, n_hidden_layers, n_nodes_per_layer = "tanh", "Adam", 4, 40
loss, metrics = "mse", ["mae", "mse"]

# Create data objects for training using scalar normalization
n_inputs = len(input_labels)
n_outputs = len(output_labels)
x = input_data
y = output_data

from idaes.core.surrogate.sampling.scaling import OffsetScaler

input_scaler = OffsetScaler.create_normalizing_scaler(x)
output_scaler = OffsetScaler.create_normalizing_scaler(y)
x = input_scaler.scale(x)
y = output_scaler.scale(y)
x = x.to_numpy()
y = y.to_numpy()

# Create data objects for training using scalar normalization
model = tf.keras.Sequential()
model.add(
    tf.keras.layers.Dense(units=n_nodes_per_layer, input_dim=len(input_labels), activation=activation))
# Create n hidden layers
for i in range(1, n_hidden_layers):
    model.add(tf.keras.layers.Dense(units=n_nodes_per_layer, activation=activation))
model.add(tf.keras.layers.Dense(units=len(output_labels)))

# Train surrogate (calls optimizer on neural network and solves for weights)
model.compile(loss=loss, optimizer=optimizer, metrics=metrics)
# Loading the weights file
mcp_save = tf.keras.callbacks.ModelCheckpoint(
    ".mdl_wts.keras", save_best_only=True, monitor="val_loss", mode="min")
history = model.fit(
    x=x, y=y, validation_split=0.2, verbose=1, epochs=2000, callbacks=[mcp_save])

from idaes.core.surrogate.keras_surrogate import KerasSurrogate
keras_surrogate = KerasSurrogate(
    model,
    input_labels=list(input_labels),
    output_labels=list(output_labels),
    input_bounds=input_bounds,
    input_scaler=input_scaler,
    output_scaler=output_scaler,
)

keras_surrogate.save_to_folder(
    keras_folder_name="temps/keras_surrogate", keras_model_name="keras_model"
)

from co2_eor.utility.sm_plotter import surrogate_parity
surrogate_parity(keras_surrogate,[training_df,validation_df],dataset_labels=["Training Data","Testing Data"],error_bound=10,filename='temps/keras_parity.pdf',show=False)
from co2_eor.utility.sm_analysis import evaluate_surrogate_metrics
metrics_df = evaluate_surrogate_metrics(keras_surrogate,[training_df,validation_df],dataset_labels=["Training Data","Testing Data"])
metrics_df.to_csv('temps/keras_metrics.csv')