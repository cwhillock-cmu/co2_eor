import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import accuracy_score, classification_report, confusion_matrix
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import (
    LabelEncoder,
    MinMaxScaler,
    RobustScaler,
    StandardScaler,
)


def learn_separating_hyperplane(
    df: pd.DataFrame,
    feature_cols: list[str],
    target_col: str,
    scaler_type: str | None = None,
    class_weight: str | None = None,
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
    # 1. Filter target columns and drop missing rows; ignore other columns
    selected_cols = feature_cols + [target_col]
    clean_df = df[selected_cols].dropna().copy()

    X = clean_df[feature_cols]
    y_raw = clean_df[target_col]

    # Encode target variable
    label_encoder = LabelEncoder()
    y = label_encoder.fit_transform(y_raw)
    class_names = [str(cls) for cls in label_encoder.classes_]

    # 2. Process Feature Scaling
    scalers = {
        "standard": StandardScaler(),
        "minmax": MinMaxScaler(),
        "robust": RobustScaler(),
    }

    scaler = None
    if scaler_type and scaler_type.lower() in scalers:
        scaler = scalers[scaler_type.lower()]
        scaler.set_output(transform="pandas")
        X_scaled_arr = scaler.fit_transform(X)
        X_proc = pd.DataFrame(X_scaled_arr, columns=feature_cols, index=X.index)
        scaling_label = scaler.__class__.__name__
    else:
        X_proc = X.copy()
        scaling_label = "None (Raw Values)"

    # 3. Fit Linear Classifier
    model = LogisticRegression(random_state=42,class_weight=class_weight)
    model.fit(X_proc, y)

    weights = model.coef_[0]
    bias = model.intercept_[0]

    # 4. Dynamic Printing of Weights, Bias, and Equations
    print("=" * 65)
    print("           MODEL WEIGHTS & HYPERPLANE PARAMETERS           ")
    print("=" * 65)
    print(f"Feature Scaling Applied : {scaling_label}")
    print(f"Bias / Intercept (b)    : {bias:+.6f}\n")

    weight_df = pd.DataFrame(
        {
            "Feature": feature_cols,
            "Weight (w)": weights,
            "Abs Weight": np.abs(weights),
        }
    ).sort_values(by="Abs Weight", ascending=False)

    print("Feature Weights Summary:")
    for _, row in weight_df.iterrows():
        print(f"  • {row['Feature']:<25} : {row['Weight (w)']:+.6f}")

    # Build primary linear equation
    terms = [f"({w:+.4f} * '{col}')" for col, w in zip(feature_cols, weights)]
    eq_str = " + ".join(terms) + f" {bias:+.4f} = 0"
    print(f"\nLearned Equation on [{scaling_label}] Space:")
    print(f"  {eq_str}\n")

    # If StandardScaler was used, derive and display the equivalent raw units equation

    if isinstance(scaler, StandardScaler):
        raw_weights = weights / scaler.scale_
        raw_bias = bias - np.sum((weights * scaler.mean_) / scaler.scale_)
        raw_terms = [
            f"({w:+.6f} * '{col}')" for col, w in zip(feature_cols, raw_weights)
        ]
        raw_eq_str = " + ".join(raw_terms) + f" {raw_bias:+.6f} = 0"
        raw_hyperplane_dict = {
            "weights":raw_weights,
            "bias":raw_bias,
        }
        print("Equivalent Raw Scale Equation (Unwound to original units):")
        print(f"  {raw_eq_str}\n")
        print(f'raw_hyperplane_dict = {raw_hyperplane_dict=}')

    # 5. Statistical Performance Analysis
    y_pred = model.predict(X_proc)
    acc = accuracy_score(y, y_pred)
    cm = confusion_matrix(y, y_pred)

    print("=" * 65)
    print("                  STATISTICAL EVALUATION                   ")
    print("=" * 65)
    print(f"Overall Accuracy Score : {acc * 100:.2f}%\n")

    cm_df = pd.DataFrame(
        cm,
        index=[f"Actual: {c}" for c in class_names],
        columns=[f"Pred: {c}" for c in class_names],
    )
    print("Confusion Matrix:")
    print(cm_df)

    print("\nDetailed Classification Report:")
    print(
        classification_report(
            y, y_pred, target_names=class_names, zero_division=0
        )
    )

    # 6. Adaptive Visualization
    n_features = len(feature_cols)

    if n_features == 1:
        fig, ax = plt.subplots(figsize=(8, 5))
        f_name = feature_cols[0]

        x_grid = np.linspace(
            X_proc[f_name].min() - 0.5, X_proc[f_name].max() + 0.5, 300
        )
        grid_df = pd.DataFrame({f_name: x_grid})
        probs = model.predict_proba(grid_df)[:, 1]

        ax.scatter(
            X_proc[f_name],
            y,
            c=y,
            cmap="coolwarm",
            edgecolors="k",
            s=60,
            zorder=3,
            label="Data Points",
        )
        ax.plot(
            x_grid,
            probs,
            color="black",
            linewidth=2,
            label="Logistic Curve $P(Y=1|X)$",
        )

        if weights[0] != 0:
            boundary_x = -bias / weights[0]
            ax.axvline(
                x=boundary_x,
                color="red",
                linestyle="--",
                linewidth=2,
                label=f"Boundary ($x = {boundary_x:.3f}$)",
            )

        ax.set_xlabel(f"{f_name} [{scaling_label}]")
        ax.set_ylabel("Probability / Class Label")
        ax.set_title(f"1D Decision Boundary for '{target_col}'")
        ax.legend()
        ax.grid(True, linestyle=":", alpha=0.6)

    elif n_features == 2:
        fig, ax = plt.subplots(figsize=(9, 6))
        f1, f2 = feature_cols

        x1_min, x1_max = X_proc[f1].min() - 0.5, X_proc[f1].max() + 0.5
        x2_min, x2_max = X_proc[f2].min() - 0.5, X_proc[f2].max() + 0.5
        xx1, xx2 = np.meshgrid(
            np.linspace(x1_min, x1_max, 250), np.linspace(x2_min, x2_max, 250)
        )

        grid_df = pd.DataFrame({f1: xx1.ravel(), f2: xx2.ravel()})
        Z = model.predict(grid_df).reshape(xx1.shape)

        ax.contourf(xx1, xx2, Z, alpha=0.25, cmap="coolwarm")
        scatter = ax.scatter(
            X_proc[f1],
            X_proc[f2],
            c=y,
            cmap="coolwarm",
            edgecolors="k",
            s=55,
            alpha=0.9,
        )

        if weights[1] != 0:
            x1_line = np.array([x1_min, x1_max])
            x2_line = -(weights[0] * x1_line + bias) / weights[1]
            ax.plot(
                x1_line,
                x2_line,
                "k--",
                linewidth=2.5,
                label="Hyperplane ($w^T x + b = 0$)",
            )

        ax.set_xlabel(f"{f1} [{scaling_label}]")
        ax.set_ylabel(f"{f2} [{scaling_label}]")
        ax.set_title(
            f"2D Decision Boundary Separation for '{target_col}' [{scaling_label}]"
        )
        ax.legend(*scatter.legend_elements(), title="Classes", loc="upper right")
        ax.grid(True, linestyle=":", alpha=0.6)

    elif n_features == 3:
        fig = plt.figure(figsize=(10, 7))
        ax = fig.add_subplot(111, projection="3d")
        f1, f2, f3 = feature_cols

        scatter = ax.scatter(
            X_proc[f1],
            X_proc[f2],
            X_proc[f3],
            c=y,
            cmap="coolwarm",
            edgecolors="k",
            s=50,
            alpha=0.85,
        )

        x1_grid = np.linspace(X_proc[f1].min(), X_proc[f1].max(), 20)
        x2_grid = np.linspace(X_proc[f2].min(), X_proc[f2].max(), 20)
        X1_mesh, X2_mesh = np.meshgrid(x1_grid, x2_grid)

        if weights[2] != 0:
            X3_mesh = -(
                weights[0] * X1_mesh + weights[1] * X2_mesh + bias
            ) / weights[2]
            ax.plot_surface(
                X1_mesh, X2_mesh, X3_mesh, alpha=0.35, cmap="coolwarm"
            )

        ax.set_xlabel(f"{f1} [{scaling_label}]")
        ax.set_ylabel(f"{f2} [{scaling_label}]")
        ax.set_zlabel(f"{f3} [{scaling_label}]")
        ax.set_title(f"3D Linear Decision Surface Plane for '{target_col}'")

    else:
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))

        scores = model.decision_function(X_proc)
        for cls_idx, cls_name in enumerate(class_names):
            ax1.hist(
                scores[y == cls_idx],
                bins=15,
                alpha=0.6,
                label=f"Class: {cls_name}",
            )
        ax1.axvline(
            0,
            color="red",
            linestyle="--",
            linewidth=2,
            label="Threshold ($w^T x + b = 0$)",
        )
        ax1.set_xlabel(f"Decision Score [{scaling_label}]")
        ax1.set_ylabel("Sample Count")
        ax1.set_title("Separation Along Hyperplane Normal Axis")
        ax1.legend()
        ax1.grid(True, linestyle=":", alpha=0.6)

        pca = PCA(n_components=2)
        X_pca = pca.fit_transform(X_proc)
        scatter = ax2.scatter(
            X_pca[:, 0],
            X_pca[:, 1],
            c=y,
            cmap="coolwarm",
            edgecolors="k",
            alpha=0.8,
        )
        ax2.set_xlabel(
            f"PCA Component 1 ({pca.explained_variance_ratio_[0]*100:.1f}%)"
        )
        ax2.set_ylabel(
            f"PCA Component 2 ({pca.explained_variance_ratio_[1]*100:.1f}%)"
        )
        ax2.set_title("2D PCA Projection of Multidimensional Space")
        ax2.grid(True, linestyle=":", alpha=0.6)

    plt.tight_layout()
    plt.show()

    pipeline = Pipeline([("scaler", scaler), ("clf", model)]) if scaler is not None else model

    return model, pipeline

# =====================================================================
# EXAMPLE USAGE DEMONSTRATION
# =====================================================================
if __name__ == "__main__":
    #df = pd.read_excel('temps/LHSdata123.xlsx')
    df = pd.read_csv('temps/LHSdata1.csv')

    # Run function specifying feature columns and 1 target column
    input_features = ['feed_co2_molflow',	
                      'feed_ch4_molflow',	
                      'feed_pressure',	
                      'feed_temperature',	
                      'recycle_solvent_molflow',
                	'HX1_outlet_temperature',	
                    'flash3_outlet_pressure',	
                    'column_stages',]
    target = "PER_ERROR"

    model,pipeline = learn_separating_hyperplane(
        df=df, feature_cols=input_features, target_col=target,scaler_type="standard",class_weight='balanced'
    )

    def filter_by_hyperplane(df, pipeline, feature_cols, drop_class):
        complete = df[feature_cols].notna().all(axis=1)
        pred = pd.Series(pipeline.predict(df.loc[complete,feature_cols]),index=df.index[complete])
        keep = complete & (pred.reindex(df.index) != drop_class)
        return df[keep].copy()

    df2 = filter_by_hyperplane(df,pipeline,input_features,1)
    df2.to_excel("temps/LHS_data1_cut.xlsx",index=False)

    model2, pipeline2 = learn_separating_hyperplane(
        df=df2, feature_cols=input_features, target_col=target,scaler_type="standard",class_weight=None
    )