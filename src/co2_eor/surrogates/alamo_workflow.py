import idaes.core
import pandas as pd
import pickle
from pathlib import Path
from contextlib import redirect_stdout

def learn_separating_hyperplane(
    df: pd.DataFrame,
    feature_cols: list[str],
    target_col: str,
    class_weight: str | None = 'balanced',
):
    """Learns a linear separating hyperplane on specified DataFrame columns with optional feature scaling.

    Parameters:
    -----------
    df : pd.DataFrame
        Source DataFrame containing features and target columns.
    feature_cols : list of str
        Column names to be used as model input features.
    target_col : str
        Column name representing the separation target label.
    scaler_type : str or None, default=None
        Type of scaling to apply to feature_cols.
        Options: 'standard' (StandardScaler), 'minmax' (MinMaxScaler),
                 'robust' (RobustScaler), or None (no scaling).

    Returns:
    --------
    model : LogisticRegression
        Fitted Scikit-Learn Logistic Regression model.
    """
    import numpy as np
    import pandas as pd
    from sklearn.decomposition import PCA
    from sklearn.linear_model import LogisticRegression
    from sklearn.metrics import accuracy_score, classification_report, confusion_matrix
    from sklearn.pipeline import Pipeline
    from sklearn.preprocessing import (
        LabelEncoder,
        StandardScaler,
    )
    # 1. Filter target columns and drop missing rows; ignore other columns
    selected_cols = feature_cols + [target_col]
    clean_df = df[selected_cols].dropna().copy()

    X = clean_df[feature_cols]
    y_raw = clean_df[target_col]

    # Encode target variable
    label_encoder = LabelEncoder()
    y = label_encoder.fit_transform(y_raw)
    class_names = [str(cls) for cls in label_encoder.classes_]

    scaler = StandardScaler()
    scaler.set_output(transform="pandas")
    X_scaled_arr = scaler.fit_transform(X)
    X_proc = pd.DataFrame(X_scaled_arr, columns=feature_cols, index=X.index)

    # 3. Fit Linear Classifier
    model = LogisticRegression(random_state=42,class_weight=class_weight)
    model.fit(X_proc, y)

    weights = model.coef_[0]
    bias = model.intercept_[0]

    raw_weights = weights / scaler.scale_
    raw_bias = bias - np.sum((weights * scaler.mean_) / scaler.scale_)
    raw_hyperplane_dict = {
        "weights":raw_weights,
        "bias":raw_bias,
    }
    print(f"  {raw_hyperplane_dict}\n")

    # 5. Statistical Performance Analysis
    y_pred = model.predict(X_proc)
    acc = accuracy_score(y, y_pred)
    cm = confusion_matrix(y, y_pred)

    cm_df = pd.DataFrame(
        cm,
        index=[f"Actual: {c}" for c in class_names],
        columns=[f"Pred: {c}" for c in class_names],
    )
    print("Confusion Matrix:")
    print(cm_df)

    pipeline = Pipeline([("scaler", scaler), ("clf", model)])

    return model, pipeline, raw_hyperplane_dict

def train_alamo_surrogate(
        full_df,
        input_labels,
        output_labels,
        input_bounds,
        save_dir:str,
        model_name:str,
        input_scales:list | None = None,
):
    input_bounds_dict = {input_labels[i]: (input_bounds[1][i],input_bounds[0][i]) for i in range(len(input_labels))}

    from idaes.core.surrogate.sampling.data_utils import split_training_validation
    training_df, validation_df = split_training_validation(full_df,0.9,seed=1)

    trainer = idaes.core.surrogate.AlamoTrainer(
        input_labels=input_labels,
        output_labels=output_labels,
        training_dataframe = training_df,
        input_bounds=input_bounds_dict,
        overwrite_files=True,
    )

    trainer.config.constant = 1
    trainer.config.linfcns = 1
    trainer.config.monomialpower = [2,3,4,5]
    trainer.config.multi2power = [1,2,3]
    trainer.config.logfcns = 1
    trainer.config.expfcns = 1
    trainer.config.ratiopower = [1,2,3]

    trainer.config.multi3power = [1,2]
    trainer.config.screener = 1
    trainer.config.modeler = 1
    trainer.config.convpen = 1
    trainer.config.maxtime = 7200
    trainer.config.ncvf = int(len(training_df))

    if input_scales is not None:
        trainer.config.xfactor = input_scales
    trainer.config.scalez = 1

    success, alm_surr, msg = trainer.train_surrogate()

    surrogate_expressions = trainer._results['Model']
    alm_surr = idaes.core.surrogate.AlamoSurrogate(surrogate_expressions, input_labels, output_labels,input_bounds_dict)
    model = alm_surr.save_to_file(save_dir+model_name+'.json', overwrite=True)

    from co2_eor.surrogates.utility.sm_plotter import surrogate_parity
    surrogate_parity(alm_surr,[training_df,validation_df],dataset_labels=["Training Data","Testing Data"],error_bound=10,filename=save_dir+model_name+'_parity.pdf',show=False)
    from co2_eor.surrogates.utility.sm_analysis import evaluate_surrogate_metrics
    metrics_df = evaluate_surrogate_metrics(alm_surr,[training_df,validation_df],dataset_labels=["Training Data","Testing Data"])
    metrics_df.to_csv(save_dir+model_name+'_metrics.csv')

if __name__ == "__main__":
    raw_df = pd.read_csv('temps/LHSdata4.csv')
    output_dir = 'temps/rectisol_surrogate/'
    Path(output_dir).mkdir(parents=True, exist_ok=True)
    co2_recovery_min = 0.9

    input_bounds = [
        [400,200,50e5,320,2500,260,10e5,10],
        [100,50,10e5,290,800,220,2e5,4]
    ]

    input_scales = [100,10,1e6,100,100,100,1e5,1]

    input_cols = ['feed_co2_molflow',	
                      'feed_ch4_molflow',	
                      'feed_pressure',	
                      'feed_temperature',	
                      'recycle_solvent_molflow',
                	'HX1_outlet_temperature',	
                    'flash3_outlet_pressure',	
                    'column_stages',]

    output_cols = ['recycle_outlet_co2_molflow',	'recycle_outlet_ch4_molflow',
                        'top_outlet_temperature',	'column_diameter',	'HX1_area',	'HX1_duty',	'pump1_volflow_in',	'pump1_work',	'comp1_work',
                                'comp2_work',	'comp2_duty',
                                        'flash1_volflow_in',	'flash2_volflow_in',	'flash3_volflow_in', 'makeup_solvent_molflow',
                                        ]
    
    cut_df = raw_df[(raw_df['recycle_outlet_co2_molflow']/raw_df['feed_co2_molflow']>co2_recovery_min) | raw_df['PER_ERROR'] == 1]

    _,_,hyperplane_dict = learn_separating_hyperplane(df=cut_df,feature_cols=input_cols,target_col='PER_ERROR',class_weight='balanced')
    hyperplane_list = [hyperplane_dict] #might want to add more cuts later
    
    clean_cut_df = cut_df[cut_df['PER_ERROR']==0]

    with open(output_dir+"alamo_output.txt", "w", encoding="utf-8") as f:
        with redirect_stdout(f):
            train_alamo_surrogate(full_df=clean_cut_df,input_labels=input_cols,output_labels=output_cols,input_bounds=input_bounds,input_scales=input_scales,save_dir=output_dir,model_name='alamo_surrogate')

    with open(output_dir+"input_bounds.pkl", "wb") as file:
        pickle.dump(input_bounds, file)

    with open(output_dir+"co2_recovery_min.pkl", "wb") as file:
        pickle.dump(co2_recovery_min, file)
 
    with open(output_dir+"hyperplane_list.pkl", "wb") as file:
        pickle.dump(hyperplane_list, file)
