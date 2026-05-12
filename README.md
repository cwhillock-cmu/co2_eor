# CO2 Enhanced Oil Recovery (EOR) Optimization Model

## Overview
This project provides a Non-Linear Programming (NLP) model built using the [IDAES](https://idaes-pse.readthedocs.io/) framework and Pyomo for optimizing CO2 Enhanced Oil Recovery (EOR) systems.

## Project Structure
* **`co2_eor/pipeline.py`**: Defines the `pipeline` unit model. It calculates pressure drop, heat transfer, and flow properties of CO2 through pipes using customized formulations (e.g., isothermal or non-isothermal heat balances).
* **`co2_eor/wellpad.py`**: Defines the `wellpad` unit model. It links CO2 injection rates to reservoir dynamics, tracking Hydrocarbon Pore Volume (HCPV), oil/gas production rates, and CO2 breakthrough using empirical curve fits and Darcy's law.
* **`co2_eor/mixer.py` & `co2_eor/splitter.py`**: Custom modified IDAES mixer and separator blocks that handle combining or splitting CO2 streams. They enforce mass, energy, and momentum balances with specific modifications like custom inequality constraints for pressures.
* **`co2_eor/flowsheets/`**: Contains the current main network assembly and optimization scripts. These scripts connect the individual unit models into full EOR distribution networks.
* **`co2_eor/util_funcs.py`**: Utility functions for defining NLP solver configurations (IPOPT, CONOPT, BONMIN, etc.) and post-processing. Includes the `export_flowsheet_to_excel` function, which automatically detects all `pipeline`, `wellpad`, and compressor blocks in a flowsheet to generate a consolidated Excel report.

*(Note: `mFTPx.py`, `thermo_config.py`, and the `deprecated` directories contain older or alternative approaches and properties.)*

## Current Workflow
The typical workflow for running an optimization, as seen in the current active flowsheets (`flowsheet5.py` and `flowsheet6.py`), follows a standard IDAES optimization pattern:

1. **Model Setup**:
   * Create a Pyomo `ConcreteModel` and an IDAES `FlowsheetBlock`.
   * Define the thermodynamic properties of CO2 using the `HelmholtzParameterBlock`.
   * Initialize costing blocks (e.g., `SSLWCosting`) for economic evaluation of the equipment.

2. **Network Assembly**:
   * Instantiate unit models (compressors, pipelines, splitters, mixers, wellpads) with their specific configuration dictionaries (e.g., length, reservoir pressure, roughness).
   * Connect the unit models using Pyomo `Arc` objects.
   * Expand the arcs to establish the full mathematical network constraints.

3. **Initialization**:
   * Fix necessary degrees of freedom (e.g., inlet mass flows, outlet pressures) to create a square problem.
   * Use `SequentialDecomposition` with a tear method (like `Wegstein`) to systematically initialize unit models in the network.
   * Propagate states and ensure mathematical feasibility before transitioning to optimization.

4. **Optimization**:
   * Define the objective function encompassing OPEX (compressor work), raw material costs, CAPEX (pipelines and compressors), and Revenue (derived from total oil production).
   * Add global constraints (e.g., a total oil production target constraint).
   * Unfix decision variables (like pipeline diameters and flow rates) and deactivate slack variables used during initialization.
   * Scale the model using Pyomo's `core.scale_model` transformation.
   * Solve the scaled Non-Linear Program using a designated solver (e.g., IPOPT or CONOPT).

5. **Post-Processing**:
   * Unscale the model and verify optimal termination.
   * Use the `export_flowsheet_to_excel` utility function from `util_funcs.py` to automatically iterate through all unit blocks (pipelines, wellpads, compressors) in the flowsheet.
   * Export the consolidated unit data and overall economics into a single Excel workbook for review.