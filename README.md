# CO2 Enhanced Oil Recovery (EOR) Optimization Model

## Overview
This project provides a Non-Linear Programming (NLP) model built using the [IDAES](https://idaes-pse.readthedocs.io/) framework and Pyomo for optimizing CO2 Enhanced Oil Recovery (EOR) systems.

## Project Structure
* **`co2_eor/MPF/`**: Contains the configuration file for building the IDAES modular properties framework Peng-Robinson Cubic equation of state. The mass state definitions didn't really work so they are not used. 
* **`co2_eor/flowsheets/`**: Contains the current main network assembly and optimization scripts. These scripts connect the individual unit models into full EOR distribution networks.
* **`co2_eor/util_funcs`**: Contains functions for exporting unit model data into dataframes and compiling them. Also has default solver options, custom initialization function (currently deprecated), and custom function for printing wide dataframes (AI-generated).

## Custom Unit Models
* **`co2_eor/gas_pipe/`**: Isothermal vapor flow pipe hydraulic model from 'Gas Pipeline Hydraulics' by Menon
* **`co2_eor/liq_pipe/`**: Nonisothermal supercritical CO2 pipe hydraulic model by Martynov et al.
* **`co2_eor/wellpattern_HH/`**: Custom well pattern model for incremental oil production from co2 injection. Currently hardcoded that inlet flow stream contains single component named 'CO2' and outlet flow streams contain 'co2' and 'ch4'. Might try to make more general later.
* **`co2_eor/processing_facility/`**: Simple yield-based separation model, also has custom costing block hardcoded into it.
* **`co2_eor/mpf_to_helmholtz/`**: Unit model that converts a MPF mixed component stream into a pure helmholtz eos stream.
* **`co2_eor/mixer_unit/`**: Standard IDAES mixer unit, no changes.
* **`co2_eor/splitter_unit/`**: Standard IDAES splitter unit, added a momentum balance type.
* **`co2_eor/heater_unit/`**: Standard IDAES heater unit, added a custom costing block that assumes heat duty will be negative.
* **`co2_eor/SSLW/`**: Standard IDAES SSLW costing file. Added smooth log functions from IDAES math util funcs. Might try to make it so these are not needed.

*(Note: `deprecated` directories contain older or alternative approaches and properties.)*

## Current Workflow
The typical workflow for running an optimization, as seen in the current active flowsheets follows a standard IDAES optimization pattern:

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
   * Solve the scaled Non-Linear Program using a designated solver (e.g., IPOPT or CONOPT). Currently CONOPT works best.

5. **Post-Processing**:
   * Unscale the model and verify optimal termination.
   * Use the `export_flowsheet_to_excel` utility function from `util_funcs.py` to automatically iterate through all unit blocks (pipelines, wellpads, compressors) in the flowsheet.
   * Export the consolidated unit data and overall economics into a single Excel workbook for review.
