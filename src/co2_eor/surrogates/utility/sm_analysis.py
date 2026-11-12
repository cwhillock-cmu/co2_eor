import numpy as np
import pandas as pd

"""
Everything in this file is AI-generated
"""

# ---------------- Default Metric Functions ---------------- #
def calculate_r2(y_true, y_pred):
    """Calculates Coefficient of Determination (R^2)."""
    ss_res = np.sum((y_true - y_pred) ** 2)
    ss_tot = np.sum((y_true - np.mean(y_true)) ** 2)
    if ss_tot == 0:
        return np.nan
    return 1.0 - (ss_res / ss_tot)


def calculate_rmse(y_true, y_pred):
    """Calculates Root Mean Squared Error (RMSE)."""
    return np.sqrt(np.mean((y_true - y_pred) ** 2))


DEFAULT_METRICS = {
    "R2": calculate_r2,
    "RMSE": calculate_rmse,
}


# ---------------- Main Evaluation Function ---------------- #
def evaluate_surrogate_metrics(
    surrogate, dataframes, dataset_labels=None, custom_metrics=None
):
    """
    Evaluates a surrogate model over a list of DataFrames and returns fit metrics
    for each output label.

    Parameters:
    -----------
    surrogate : object
        Surrogate object with `input_labels()`, `output_labels()`, and `evaluate_surrogate()`.
    dataframes : list of pd.DataFrame
        List of dataframes to evaluate.
    dataset_labels : list of str, optional
        Column labels for each dataframe. Defaults to ["data 1", "data 2", ...] if None.
    custom_metrics : dict, optional
        Dictionary mapping metric names (str) to functions with signature `fn(y_true, y_pred)`.
        These will be merged with or override DEFAULT_METRICS.

    Returns:
    --------
    pd.DataFrame
        Pandas DataFrame where columns are `dataset_labels` and rows are indexed
        by (Output Label, Metric).
    """
    # ---------------- Input Validation ---------------- #
    if not isinstance(dataframes, (list, tuple)) or len(dataframes) == 0:
        raise ValueError("'dataframes' must be a non-empty list or tuple of pandas DataFrames.")

    for idx, df in enumerate(dataframes):
        if not isinstance(df, pd.DataFrame):
            raise TypeError(f"Element at index {idx} in 'dataframes' is not a pandas DataFrame.")

    if dataset_labels is None:
        dataset_labels = [f"data {i + 1}" for i in range(len(dataframes))]
    else:
        if not isinstance(dataset_labels, (list, tuple)):
            raise TypeError("'dataset_labels' must be a list or tuple of strings.")
        if len(dataset_labels) != len(dataframes):
            raise ValueError(
                f"Length of 'dataset_labels' ({len(dataset_labels)}) does not match "
                f"length of 'dataframes' ({len(dataframes)})."
            )

    # ---------------- Metric Modularity Setup ---------------- #
    metrics_to_compute = DEFAULT_METRICS.copy()
    if custom_metrics is not None:
        if not isinstance(custom_metrics, dict):
            raise TypeError("'custom_metrics' must be a dictionary mapping metric names to functions.")
        metrics_to_compute.update(custom_metrics)

    input_labels = surrogate.input_labels()
    output_labels = surrogate.output_labels()

    # Validate dataframe columns
    for idx, df in enumerate(dataframes):
        missing_inputs = [col for col in input_labels if col not in df.columns]
        missing_outputs = [col for col in output_labels if col not in df.columns]
        if missing_inputs or missing_outputs:
            raise KeyError(
                f"DataFrame at index {idx} is missing required columns. "
                f"Missing inputs: {missing_inputs}, Missing outputs: {missing_outputs}"
            )

    # ---------------- Metric Calculations ---------------- #
    results = {}

    for label, df in zip(dataset_labels, dataframes):
        input_data = df[input_labels]
        output_true = df[output_labels]
        output_pred = surrogate.evaluate_surrogate(input_data)

        df_metrics = {}
        for out_col in output_labels:
            y_true = output_true[out_col].values
            y_pred = output_pred[out_col].values

            for metric_name, metric_fn in metrics_to_compute.items():
                df_metrics[(out_col, metric_name)] = metric_fn(y_true, y_pred)

        results[label] = df_metrics

    # Construct formatted DataFrame with dataset_labels as columns
    metrics_df = pd.DataFrame(results)
    metrics_df.index = pd.MultiIndex.from_tuples(
        metrics_df.index, names=["Output Label", "Metric"]
    )

    return metrics_df