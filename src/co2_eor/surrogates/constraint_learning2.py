import cvxpy as cp
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA
from sklearn.metrics import accuracy_score, classification_report, confusion_matrix
from sklearn.preprocessing import (
    LabelEncoder,
    MinMaxScaler,
    RobustScaler,
    StandardScaler,
)


def learn_polyhedral_classifier_cvxpy(
    df: pd.DataFrame,
    feature_cols: list[str],
    target_col: str,
    k_inequalities: int = 4,
    scaler_type: str | None = None,
    max_iters: int = 5,
    C_pos: float = 1.0,
    C_neg: float = 1.0,
):
    """Learns a system of K simultaneous linear inequalities (Ax <= b) using CVXPY.

    Parameters:
    -----------
    df : pd.DataFrame
        Source DataFrame containing features and target column.
    feature_cols : list of str
        Column names to be used as model input features.
    target_col : str
        Column name representing the binary separation target.
    k_inequalities : int, default=4
        Number of simultaneous linear inequalities (hyperplanes) to learn.
    scaler_type : str or None, default=None
        Scaling method ('standard', 'minmax', 'robust', or None).
    max_iters : int, default=5
        Maximum alternating optimization iterations for CVXPY.
    C_pos, C_neg : float, default=1.0
        Soft-margin penalty parameters for positive and negative class slacks.
    """
    # 1. Isolate target columns and handle missing values
    selected_cols = feature_cols + [target_col]
    clean_df = df[selected_cols].dropna().copy()

    X_raw = clean_df[feature_cols].values
    y_raw = clean_df[target_col].values

    # Encode target labels into binary (1 = inside polyhedron, 0 = outside)
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
        X_scaled = scaler.fit_transform(X_raw)
        scaling_label = scaler.__class__.__name__
    else:
        X_scaled = X_raw.copy()
        scaling_label = "None (Raw Values)"

    # Separate positive (Class 1) and negative (Class 0) data points
    pos_mask = y == 1
    neg_mask = y == 0

    X_pos = X_scaled[pos_mask]
    X_neg = X_scaled[neg_mask]

    N_pos, d = X_pos.shape
    N_neg = X_neg.shape[0]

    # 3. Initialize Hyperplane Directions using K-Means on Negative Data
    c_pos = X_pos.mean(axis=0)

    if N_neg >= k_inequalities:
        kmeans = KMeans(
            n_clusters=k_inequalities, random_state=42, n_init=10
        ).fit(X_neg)
        dirs = kmeans.cluster_centers_ - c_pos
        norm_dirs = np.linalg.norm(dirs, axis=1, keepdims=True)
        norm_dirs[norm_dirs == 0] = 1.0
        dirs = dirs / norm_dirs
    else:
        dirs = np.random.randn(k_inequalities, d)
        dirs /= np.linalg.norm(dirs, axis=1, keepdims=True)

    A_init = dirs
    b_init = np.sum(A_init * c_pos, axis=1) + 1.0

    # Initial assignment of negative points to nearest boundary face
    scores = X_neg @ A_init.T - b_init
    assignments = np.argmax(scores, axis=1)

    A_learned = A_init
    b_learned = b_init

    # 4. Alternating Optimization with CVXPY
    for iteration in range(max_iters):
        A_var = cp.Variable((k_inequalities, d))
        b_var = cp.Variable(k_inequalities)
        xi = cp.Variable((N_pos, k_inequalities), nonneg=True)
        eta = cp.Variable(N_neg, nonneg=True)

        constraints = []

        # Enforce all positive points to lie INSIDE all K inequalities
        for k in range(k_inequalities):
            constraints.append(
                X_pos @ A_var[k, :] - b_var[k] <= -1.0 + xi[:, k]
            )

        # Enforce negative points to lie OUTSIDE their assigned inequality face
        for k in range(k_inequalities):
            idx_k = np.where(assignments == k)[0]
            if len(idx_k) > 0:
                constraints.append(
                    X_neg[idx_k, :] @ A_var[k, :] - b_var[k]
                    >= 1.0 - eta[idx_k]
                )

        # Soft-margin convex objective function
        objective = cp.Minimize(
            0.5 * cp.sum_squares(A_var)
            + C_pos * cp.sum(xi)
            + C_neg * cp.sum(eta)
        )

        problem = cp.Problem(objective, constraints)

        # Solve QP via CVXPY solvers
        try:
            problem.solve(solver=cp.ECOS,verbose=False)
        except Exception:
            problem.solve()

        if A_var.value is None or b_var.value is None:
            break

        A_learned = A_var.value
        b_learned = b_var.value

        # Re-assign negative points to the face they violate the most
        scores = X_neg @ A_learned.T - b_learned
        new_assignments = np.argmax(scores, axis=1)
        if np.array_equal(new_assignments, assignments):
            break
        assignments = new_assignments

    # 5. Print Weights, Biases, and Explicit Linear Inequalities
    print("=" * 70)
    print("      CVXPY POLYHEDRAL CLASSIFIER (SIMULTANEOUS INEQUALITIES)     ")
    print("=" * 70)
    print(f"Feature Scaling Applied   : {scaling_label}")
    print(f"Number of Inequalities (K): {k_inequalities}\n")

    print(f"Learned Polyhedral Inequalities on [{scaling_label}] Space:")
    print("  Condition: Point x is in Class 1 if ALL inequalities hold (Ax <= b)\n")

    for k in range(k_inequalities):
        terms = [
            f"({A_learned[k, j]:+.4f} * '{feature_cols[j]}')"
            for j in range(d)
        ]
        eq_str = " + ".join(terms) + f" <= {b_learned[k]:+.4f}"
        print(f"  [Inequality {k+1}]:  {eq_str}")

    # Unwind equations if StandardScaler was applied
    if isinstance(scaler, StandardScaler):
        A_raw = A_learned / scaler.scale_
        b_raw = b_learned + np.sum(
            (A_learned * scaler.mean_) / scaler.scale_, axis=1
        )
        print("\nEquivalent Raw Units Inequalities (Unwound to Original Scale):")
        for k in range(k_inequalities):
            raw_terms = [
                f"({A_raw[k, j]:+.6f} * '{feature_cols[j]}')" for j in range(d)
            ]
            raw_eq_str = " + ".join(raw_terms) + f" <= {b_raw[k]:+.6f}"
            print(f"  [Inequality {k+1}]:  {raw_eq_str}")

    # 6. Evaluate Predictions & Statistical Performance
    # Point belongs to Class 1 if max_k (A_k x - b_k) <= 0
    max_scores = np.max(X_scaled @ A_learned.T - b_learned, axis=1)
    y_pred = (max_scores <= 0).astype(int)

    acc = accuracy_score(y, y_pred)
    cm = confusion_matrix(y, y_pred)

    print("\n" + "=" * 70)
    print("                  STATISTICAL EVALUATION                   ")
    print("=" * 70)
    print(f"Overall Classification Accuracy: {acc * 100:.2f}%\n")

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

    # 7. Adaptive Visualization
    if d == 1:
        # 1D Bounding Interval
        fig, ax = plt.subplots(figsize=(8, 5))
        f1 = feature_cols[0]

        ax.scatter(
            X_scaled[:, 0],
            y,
            c=y,
            cmap="coolwarm",
            edgecolors="k",
            s=60,
            zorder=3,
        )

        for k in range(k_inequalities):
            if A_learned[k, 0] != 0:
                bound_x = b_learned[k] / A_learned[k, 0]
                ax.axvline(
                    x=bound_x,
                    color="black",
                    linestyle="--",
                    label=f"Bound {k+1} ($x = {bound_x:.2f}$)",
                )

        ax.set_xlabel(f"{f1} [{scaling_label}]")
        ax.set_ylabel("Class Label")
        ax.set_title(
            f"1D Bounding Inequalities ({k_inequalities} Thresholds) for '{target_col}'"
        )
        ax.legend()
        ax.grid(True, linestyle=":", alpha=0.6)

    elif d == 2:
        # 2D Polyhedral Region & Intersecting Lines
        fig, ax = plt.subplots(figsize=(9, 6))
        f1, f2 = feature_cols

        x1_min, x1_max = (
            X_scaled[:, 0].min() - 0.5,
            X_scaled[:, 0].max() + 0.5,
        )
        x2_min, x2_max = (
            X_scaled[:, 1].min() - 0.5,
            X_scaled[:, 1].max() + 0.5,
        )
        xx1, xx2 = np.meshgrid(
            np.linspace(x1_min, x1_max, 300), np.linspace(x2_min, x2_max, 300)
        )

        grid_pts = np.c_[xx1.ravel(), xx2.ravel()]
        grid_max_scores = np.max(
            grid_pts @ A_learned.T - b_learned, axis=1
        ).reshape(xx1.shape)
        grid_pred = (grid_max_scores <= 0).astype(int)

        ax.contourf(xx1, xx2, grid_pred, alpha=0.25, cmap="coolwarm")
        scatter = ax.scatter(
            X_scaled[:, 0],
            X_scaled[:, 1],
            c=y,
            cmap="coolwarm",
            edgecolors="k",
            s=50,
            alpha=0.9,
        )

        # Plot individual linear boundary lines
        x1_vals = np.array([x1_min, x1_max])
        for k in range(k_inequalities):
            if A_learned[k, 1] != 0:
                x2_vals = (b_learned[k] - A_learned[k, 0] * x1_vals) / A_learned[
                    k, 1
                ]
                ax.plot(
                    x1_vals,
                    x2_vals,
                    "--",
                    linewidth=1.8,
                    label=f"Line {k+1}: $A_{{{k+1}}} x = b_{{{k+1}}}$",
                )

        ax.set_xlim(x1_min, x1_max)
        ax.set_ylim(x2_min, x2_max)
        ax.set_xlabel(f"{f1} [{scaling_label}]")
        ax.set_ylabel(f"{f2} [{scaling_label}]")
        ax.set_title(
            f"2D Bounded Polyhedron ({k_inequalities} Intersecting Inequalities)"
        )
        ax.legend(loc="upper right")
        ax.grid(True, linestyle=":", alpha=0.6)

    elif d == 3:
        # 3D Polyhedral Planes
        fig = plt.figure(figsize=(10, 7))
        ax = fig.add_subplot(111, projection="3d")
        f1, f2, f3 = feature_cols

        ax.scatter(
            X_scaled[:, 0],
            X_scaled[:, 1],
            X_scaled[:, 2],
            c=y,
            cmap="coolwarm",
            edgecolors="k",
            s=45,
            alpha=0.8,
        )

        x1_g = np.linspace(X_scaled[:, 0].min(), X_scaled[:, 0].max(), 15)
        x2_g = np.linspace(X_scaled[:, 1].min(), X_scaled[:, 1].max(), 15)
        X1_m, X2_m = np.meshgrid(x1_g, x2_g)

        for k in range(k_inequalities):
            if A_learned[k, 2] != 0:
                X3_m = (
                    b_learned[k]
                    - A_learned[k, 0] * X1_m
                    - A_learned[k, 1] * X2_m
                ) / A_learned[k, 2]
                ax.plot_surface(X1_m, X2_m, X3_m, alpha=0.2)

        ax.set_xlabel(f1)
        ax.set_ylabel(f2)
        ax.set_zlabel(f3)
        ax.set_title(
            f"3D Polyhedral Bounding Planes ({k_inequalities} Inequalities)"
        )

    else:
        # >3D Feature Space Visualization
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))

        # Max constraint score distribution
        for cls_idx, cls_name in enumerate(class_names):
            ax1.hist(
                max_scores[y == cls_idx],
                bins=15,
                alpha=0.6,
                label=f"Class: {cls_name}",
            )
        ax1.axvline(
            0,
            color="red",
            linestyle="--",
            linewidth=2,
            label=r"Polyhedron Bound ($\max (Ax - b) = 0$)",
        )
        ax1.set_xlabel(r"Maximum Inequality Value $\max_k (A_k x - b_k)$")
        ax1.set_ylabel("Sample Count")
        ax1.set_title("Margin Distribution Relative to Polyhedron")
        ax1.legend()
        ax1.grid(True, linestyle=":", alpha=0.6)

        # 2D PCA projection
        pca = PCA(n_components=2)
        X_pca = pca.fit_transform(X_scaled)
        ax2.scatter(
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
        ax2.set_title("2D PCA Projection of Polyhedral Space")
        ax2.grid(True, linestyle=":", alpha=0.6)

    plt.tight_layout()
    plt.show()

    return A_learned, b_learned


# =====================================================================
# EXAMPLE DEMONSTRATION
# =====================================================================
if __name__ == "__main__":
    df = pd.read_excel('temps/LHSdata123.xlsx')

    # Run function specifying feature columns and 1 target column
    input_features = ['feed_co2_molflow',	'feed_ch4_molflow',	'feed_pressure',	'feed_temperature',	'recycle_solvent_molflow',
                	'HX1_outlet_temperature',	'flash3_outlet_pressure',	'column_stages',]
    target = "PER_ERROR"

    # 2. Run CVXPY Polyhedral Classifier learning simultaneous linear inequalities
    A_mat, b_vec = learn_polyhedral_classifier_cvxpy(
        df=df,
        feature_cols=input_features,
        target_col=target,
        k_inequalities=4,
        scaler_type="standard",
        max_iters=50,
    )