#################################################################################
# The Institute for the Design of Advanced Energy Systems Integrated Platform
# Framework (IDAES IP) was produced under the DOE Institute for the
# Design of Advanced Energy Systems (IDAES).
#
# Copyright (c) 2018-2026 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory,
# National Technology & Engineering Solutions of Sandia, LLC, Carnegie Mellon
# University, West Virginia University Research Corporation, et al.
# All rights reserved.  Please see the files COPYRIGHT.md and LICENSE.md
# for full copyright and license information.
#################################################################################
"""
Standard plotting methods (scatter, parity and residual plots).

Author: B. Paul
Modified for co2_eor
"""

# TODO: Missing docstrings
# pylint: disable=missing-function-docstring

# Import statements
from itertools import combinations  # used to pick xvar pairs in scatter3D

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

# TODO: Consider adding seaborn in future
from matplotlib.backends.backend_pdf import PdfPages


def surrogate_scatter2D(surrogate, dataframe, filename=None, show=True):
    input_data = dataframe[surrogate.input_labels()]
    output_data = dataframe[surrogate.output_labels()]
    output_surrogate = surrogate.evaluate_surrogate(input_data)
    return _scatter2D(
        xdata=input_data.values,
        zdata=output_data.values,
        zfit=output_surrogate.values,
        xlabels=surrogate.input_labels(),
        zlabels=surrogate.output_labels(),
        show=show,
        filename=filename,
    )


def _scatter2D(
    xdata, zdata, zfit, xlabels=None, zlabels=None, show=True, filename=None
):
    """
    Plots training data as 2D scatter, with overlaid model output data.

    Required: input variable data, output variable data, model fit data

    Optional: input and output labels (will be generated as needed), filename
    to save plots to PDF

    Generally, x (input) variables are indexed by 'i', and z (output) variables
    are indexed by 'j'.
    """

    numins = np.shape(xdata)[1]  # number of input variables
    numouts = np.shape(zdata)[1]  # number of output variables

    # check and add labels if missing
    if xlabels is None:
        xlabels = []
        for i in range(numins):
            xlabels.append("x" + str(i + 1))
    if zlabels is None:
        zlabels = []
        for j in range(numouts):
            zlabels.append("z" + str(j + 1))

    if filename is not None:
        if filename[-4:] != ".pdf":  # checking if extension is present
            # PYLINT-TODO
            # pylint: disable-next=broad-exception-raised
            raise Exception(
                "Filename does not end in .pdf, please amend "
                "filename with the proper extension."
            )
        pdf = PdfPages(filename)

    fig_list = []

    for j in range(numouts):  # loop over all outputs, zj
        for i in range(numins):  # plot every possible zj = f(xi)
            fig = plt.figure()
            ax = fig.add_subplot()
            ax.scatter(xdata[:, i], zdata[:, j], c="grey", marker="s", label="Data")
            ax.scatter(xdata[:, i], zfit[:, j], c="b", marker=".", label="Model")

            ax.set_xlabel(xlabels[i])
            ax.set_ylabel(zlabels[j])
            ax.set_title("2D Scatter Plot")
            ax.legend()

            if show is True:
                plt.show()
            if filename is not None:
                pdf.savefig(fig)

            fig_list.append(fig)

    if filename is not None:  # place outside loop to avoid closing/reopening
        pdf.close()

    return fig_list


def surrogate_scatter3D(surrogate, dataframe, filename=None, show=True):
    input_data = dataframe[surrogate.input_labels()]
    output_data = dataframe[surrogate.output_labels()]
    output_surrogate = surrogate.evaluate_surrogate(input_data)
    return _scatter3D(
        xdata=input_data.values,
        zdata=output_data.values,
        zfit=output_surrogate.values,
        xlabels=surrogate.input_labels(),
        zlabels=surrogate.output_labels(),
        show=show,
        filename=filename,
    )


def _scatter3D(
    xdata, zdata, zfit, xlabels=None, zlabels=None, show=True, filename=None
):
    """
    Plots training data as 3D scatter, with overlaid model output data.

    Required: input variable data, output variable data, model fit data

    Optional: input and output labels (will be generated as needed), filename
    to save plots to PDF

    Generally, x (input) variables are indexed by 'i', and z (output) variables
    are indexed by 'j'.
    """

    numins = np.shape(xdata)[1]  # number of input variables
    numouts = np.shape(zdata)[1]  # number of output variables

    # check and add labels if missing
    if xlabels is None:
        xlabels = []
        for i in range(numins):
            xlabels.append("x" + str(i + 1))
    if zlabels is None:
        zlabels = []
        for j in range(numouts):
            zlabels.append("z" + str(j + 1))

    if filename is not None:
        if filename[-4:] != ".pdf":  # checking if extension is present
            # PYLINT-TODO
            # pylint: disable-next=broad-exception-raised
            raise Exception(
                "Filename does not end in .pdf, please amend "
                "filename with the proper extension."
            )
        pdf = PdfPages(filename)

    fig_list = []

    for j in range(numouts):  # loop over all outputs, zj
        for pair in list(combinations(range(numins), 2)):  # pick two x vars
            a, b = pair[0], pair[1]  # indices for the x variables picked
            fig = plt.figure()
            ax = fig.add_subplot(projection="3d")

            ax.scatter(
                xdata[:, a],
                xdata[:, b],
                zdata[:, j],
                c="grey",
                marker="s",
                label="Data",
            )
            ax.scatter(
                xdata[:, a], xdata[:, b], zfit[:, j], c="b", marker=".", label="Model"
            )
            ax.set_xlabel(xlabels[a])
            ax.set_ylabel(xlabels[b])
            ax.set_zlabel(zlabels[j])
            ax.set_title("3D Scatter Plot")
            ax.legend()

            if show is True:
                plt.show()
            if filename is not None:
                pdf.savefig(fig)

            fig_list.append(fig)

    if filename is not None:  # place outside loop to avoid closing/reopening
        pdf.close()

    return fig_list

def surrogate_parity(
    surrogate,
    dataframes,
    dataset_labels=None,
    filename=None,
    show=True,
    colors=None,
    error_bound=0,
):
    """
    THIS FUNCTION IS AI GENERATED
    Evaluates a surrogate model over a list of DataFrames and generates multi-dataset
    parity plots with custom or default dataset labels, colors, and optional percentage error bound lines.

    Parameters:
    -----------
    surrogate : object
        Surrogate object with methods `input_labels()`, `output_labels()`, and `evaluate_surrogate()`.
    dataframes : list of pd.DataFrame
        List of dataframes containing input and target output variables.
    dataset_labels : list of str, optional
        Legend labels corresponding to each dataframe (e.g., ["training points", "testing points"]).
        Defaults to ["data 1", "data 2", ...] if omitted.
    filename : str, optional
        PDF filename to save the generated plots to.
    show : bool, optional
        Whether to call plt.show() for each figure (default is True).
    colors : list of str, optional
        Custom list of colors for each dataframe's scatter points.
    error_bound : float or int, optional
        Percentage error bound x (e.g., passing 10 plots ±10% dashed lines).
        Passing 0 or omitting it (default 0) disables drawing the error bound lines.
    """
    # ------------------ Input Validation ------------------ #
    if not isinstance(dataframes, (list, tuple)) or len(dataframes) == 0:
        raise ValueError("'dataframes' must be a non-empty list or tuple of pandas DataFrames.")

    for idx, df in enumerate(dataframes):
        if not isinstance(df, pd.DataFrame):
            raise TypeError(f"Element at index {idx} in 'dataframes' is not a pandas DataFrame.")

    # Validate dataset_labels
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

        for idx, label in enumerate(dataset_labels):
            if not isinstance(label, str):
                raise TypeError(f"Label at index {idx} in 'dataset_labels' must be a string.")

    # Validate error_bound parameter
    if error_bound is not None:
        if isinstance(error_bound, bool) or not isinstance(error_bound, (int, float)):
            raise TypeError("'error_bound' must be a numeric value (int or float).")
        if error_bound < 0:
            raise ValueError("'error_bound' must be a non-negative number.")

    input_labels = surrogate.input_labels()
    output_labels = surrogate.output_labels()

    for idx, df in enumerate(dataframes):
        missing_inputs = [col for col in input_labels if col not in df.columns]
        missing_outputs = [col for col in output_labels if col not in df.columns]
        if missing_inputs or missing_outputs:
            raise KeyError(
                f"DataFrame at index {idx} is missing required columns. "
                f"Missing inputs: {missing_inputs}, Missing outputs: {missing_outputs}"
            )

    if colors is not None:
        if not isinstance(colors, (list, tuple)) or len(colors) != len(dataframes):
            raise ValueError("Length of 'colors' list must match length of 'dataframes'.")

    # ------------------ Data Processing ------------------ #
    zdata_list = []
    zfit_list = []

    for df in dataframes:
        input_data = df[input_labels]
        output_data = df[output_labels]
        output_surrogate = surrogate.evaluate_surrogate(input_data)

        zdata_list.append(output_data.values)
        zfit_list.append(output_surrogate.values)

    return _parity(
        zdata_list=zdata_list,
        zfit_list=zfit_list,
        dataset_labels=dataset_labels,
        zlabels=output_labels,
        colors=colors,
        error_bound=error_bound,
        show=show,
        filename=filename,
    )


def _parity(
    zdata_list,
    zfit_list,
    dataset_labels,
    zlabels=None,
    colors=None,
    error_bound=0,
    show=True,
    filename=None,
):
    """
    Plots model output against data output for multiple datasets on the same figures.
    """
    numouts = np.shape(zdata_list[0])[1]

    if zlabels is None:
        zlabels = [f"z{j + 1}" for j in range(numouts)]

    if filename is not None:
        if filename[-4:] != ".pdf":
            raise ValueError(
                "Filename does not end in .pdf, please amend "
                "filename with the proper extension."
            )
        pdf = PdfPages(filename)

    if colors is None:
        default_palette = ["blue", "red", "green", "purple", "orange", "cyan", "magenta"]
        colors = [default_palette[i % len(default_palette)] for i in range(len(zdata_list))]

    fig_list = []

    for j in range(numouts):
        all_mins = []
        all_maxs = []
        for zdata, zfit in zip(zdata_list, zfit_list):
            all_mins.extend([np.min(zdata[:, j]), np.min(zfit[:, j])])
            all_maxs.extend([np.max(zdata[:, j]), np.max(zfit[:, j])])

        minzj = min(all_mins)
        maxzj = max(all_maxs)

        fig = plt.figure()
        ax = fig.add_subplot()

        # Solid 1:1 Parity line
        ax.plot([minzj, maxzj], [minzj, maxzj], ls="-", c="grey", label="Parity")

        # Optional dashed +-x% error bound lines
        if error_bound and error_bound > 0:
            y_upper = [minzj * (1 + error_bound / 100.0), maxzj * (1 + error_bound / 100.0)]
            y_lower = [minzj * (1 - error_bound / 100.0), maxzj * (1 - error_bound / 100.0)]
            ax.plot([minzj, maxzj], y_upper, ls="--", c="grey", label=f"±{error_bound}% Bounds")
            ax.plot([minzj, maxzj], y_lower, ls="--", c="grey")

        # Plot dataset points
        for k, (zdata, zfit) in enumerate(zip(zdata_list, zfit_list)):
            ax.scatter(
                zdata[:, j],
                zfit[:, j],
                s=12,
                c=colors[k],
                marker=".",
                label=dataset_labels[k],
            )

        ax.set_xlim(minzj, maxzj)
        ax.set_ylim(minzj, maxzj)
        ax.set_xlabel(zlabels[j])
        ax.set_ylabel("Model Output")
        ax.set_title(f"Parity Plot - {zlabels[j]}")
        ax.legend()

        if show:
            plt.show()
        if filename is not None:
            pdf.savefig(fig)

        fig_list.append(fig)

    if filename is not None:
        pdf.close()

    return fig_list

def surrogate_residual(
    surrogate, dataframe, filename=None, relative_error=False, show=True
):
    input_data = dataframe[surrogate.input_labels()]
    output_data = dataframe[surrogate.output_labels()]
    output_surrogate = surrogate.evaluate_surrogate(input_data)
    residual = np.abs(output_data - output_surrogate)
    if relative_error is True:
        residual = np.divide(residual, np.maximum(output_data, 1.0))
    return _residual(
        xdata=input_data.values,
        residual=residual.values,
        xlabels=surrogate.input_labels(),
        elabels=surrogate.output_labels(),
        show=show,
        filename=filename,
    )


def _residual(xdata, residual, xlabels=None, elabels=None, show=True, filename=None):
    """
    Plots model residual against data output.

    Required: input variable data, model error data

    Optional: input and error labels (will be generated as needed), filename
    to save plots to PDF

    Generally, x (input) variables are indexed by 'i', and z (output) variables
    are indexed by 'j'.
    """
    numins = np.shape(xdata)[1]  # number of input variables
    numouts = np.shape(residual)[1]  # number of output variables

    # check and add labels if missing

    if xlabels is None:
        xlabels = []
        for j in range(numins):
            xlabels.append("x" + str(j + 1))
    if elabels is None:
        for i in range(numouts):
            elabels.append("z" + str(j + 1) + " Error")

    if filename is not None:
        if filename[-4:] != ".pdf":  # checking if extension is present
            # PYLINT-TODO
            # pylint: disable-next=broad-exception-raised
            raise Exception(
                "Filename does not end in .pdf, please amend "
                "filename with the proper extension."
            )
        pdf = PdfPages(filename)

    fig_list = []

    for i in range(numins):
        for j in range(numouts):  # loop over all outputs, zj
            fig = plt.figure()
            ax = fig.add_subplot()
            ax.scatter(xdata[:, i], residual[:, j], marker=".", label=elabels[j])
            ax.set_xlabel(xlabels[i])
            ax.set_ylabel(elabels[j])
            ax.set_title("Residual Plot")
            ax.legend()

            if show is True:
                plt.show()
            if filename is not None:
                pdf.savefig(fig)

            fig_list.append(fig)

    if filename is not None:  # place outside loop to avoid closing/reopening
        pdf.close()

    return fig_list
