from pyomo.environ import units as units
from idaes.core import Component, PhaseType, LiquidPhase, VaporPhase
from idaes.models.properties.modular_properties.state_definitions import FTPx
from idaes.models.properties.modular_properties.eos.ceos import Cubic, CubicType
from idaes.models.properties.modular_properties.pure import NIST
from idaes.models.properties.modular_properties.pure import ChungViscosityPure
from idaes.models.properties.modular_properties.pure.ChapmanEnskog import collision_integral_neufeld_callback
from idaes.models.properties.modular_properties.transport_properties import ViscosityWilke, NoMethod
from idaes.models.properties.modular_properties.transport_properties.viscosity_wilke import wilke_phi_ij_callback, herring_zimmer_phi_ij_callback
from co2_eor import mFTPx, mFPhx

configuration = {
    "base_units":{
        "time":units.s,
        "length":units.m,
        "mass":units.kg,
        "amount":units.mol,
        "temperature":units.K,
    },
    "components":{
        "co2":{
            "type":Component,
            "elemental_composition":{"C":1,"O":2},
            "enth_mol_ig_comp":NIST,
            "entr_mol_ig_comp":NIST,
            "cp_mol_ig_comp":NIST,
            "pressure_sat_comp":NIST,
            "valid_phase_types":[PhaseType.vaporPhase],
            "visc_d_phase_comp":{"Vap":ChungViscosityPure},
            #"viscosity_collision_integral_callback":collision_integral_neufeld_callback,
            "parameter_data":{
                "mw":(44.009e-3,units.kg/units.mol),
                "pressure_crit":(73.825*100000,units.Pa),
                "temperature_crit":(304.23,units.K),
                "omega":0.225,
                "compress_fact_crit":0.275,
                "dens_mol_crit":(10590,units.mol/units.m**3),
                "cp_mol_ig_comp_coeff": {
                    "A":(24.99735,units.J/units.mol/units.K),
                    "B":(55.18696,units.J/units.mol/units.K**2),
                    "C":(-33.69137,units.J/units.mol/units.K**3),
                    "D":(7.948387,units.J/units.mol/units.K**4),
                    "E":(-0.136638,units.J*units.K**2/units.mol/units.K),
                    "F":(-403.6075,units.kJ/units.mol),
                    "G":(228.2431,units.J/units.mol/units.K),
                    "H":(-393.5224,units.kJ/units.mol),
                },
                "dipole_moment":(0,units.debye),
                "association_factor_chung":0,
                "pressure_sat_comp_coeff":{
                        "A":6.81228,
                        "B":(1301.679,units.K),
                        "C":(-3.494,units.K),
                    },
                "entr_mol_form_liq_comp_ref":(0,units.J/units.mol/units.K), #double check this
                "enth_mol_form_liq_comp_ref":(0,units.J/units.mol), #double check this
            },
        },
        "ch4":{
            "type":Component,
            "elemental_composition":{"C":1,"H":4},
            "enth_mol_ig_comp":NIST,
            "entr_mol_ig_comp":NIST,
            "cp_mol_ig_comp":NIST,
            "pressure_sat_comp":NIST,
            "valid_phase_types":[PhaseType.vaporPhase],
            "visc_d_phase_comp":{"Vap":ChungViscosityPure},
            #"viscosity_collision_integral_callback":collision_integral_neufeld_callback,
            "parameter_data":{
                "mw":(16.0425e-3,units.kg/units.mol),
                "pressure_crit":(46.1*100000,units.Pa),
                "temperature_crit":(190.6,units.K),
                "omega":0.01142, #https://coolprop.org/fluid_properties/fluids/Methane.html
                "dens_mol_crit":(10139,units.mol/units.m**3),
                "compress_fact_crit":0.2869, #P/rho/R/T
                "cp_mol_ig_comp_coeff": {
                    "A":(-0.703029,units.J/units.mol/units.K),
                    "B":(108.4773,units.J/units.mol/units.K**2),
                    "C":(-42.52157,units.J/units.mol/units.K**3),
                    "D":(5.862788,units.J/units.mol/units.K**4),
                    "E":(0.678565,units.J*units.K**2/units.mol/units.K),
                    "F":(-76.84376,units.kJ/units.mol),
                    "G":(158.7163,units.J/units.mol/units.K),
                    "H":(-74.87310,units.kJ/units.mol),
                },
                "dipole_moment":(0,units.debye),
                "association_factor_chung":0,
                "pressure_sat_comp_coeff":{ #really shouldn't need these
                        "A":4.22061,
                        "B":(516.689,units.K),
                        "C":(11.223,units.K),
                    },
                "entr_mol_form_liq_comp_ref":(0,units.J/units.mol/units.K), #double check this
                "enth_mol_form_liq_comp_ref":(0,units.J/units.mol), #double check this
            },
        },
    },
    "phases":{
        #"Liq":{
        #    "type":LiquidPhase,
        #    "equation_of_state":Cubic,
        #    "equation_of_state_options":{"type":CubicType.PR},
        #    "visc_d_phase": NoMethod,
        #},
        "Vap":{
            "type":VaporPhase,
            "equation_of_state":Cubic,
            "equation_of_state_options":{"type":CubicType.PR},
            "visc_d_phase": ViscosityWilke,
            "transport_property_options": {"viscosity_phi_ij_callback": wilke_phi_ij_callback,}
        },
    },
    "state_definition":mFTPx,
    "state_bounds":{
        "flow_mass":(-20000,1,20000,units.kg/units.s),
        "temperature":(216,298,1500,units.K),
        "pressure":(5e4,80*100000,2000*100000,units.Pa)
    },
    "pressure_ref":(101325,units.Pa), #double check this
    "temperature_ref":(298.15,units.K), #double check this
    "parameter_data":{ #https://doi.org/10.1021/ie301012q 
        "PR_kappa":{
            ("co2","co2"):0,
            ("co2","ch4"):0.1094, #TRIPLE CHECK THESE NUMBERS!!!
            ("ch4","ch4"):0,
            ("ch4","co2"):0.1094,
        }
    }
}

"""
    "state_definition":mFTPx,
    "state_bounds":{
        "flow_mass":(-20000,1,20000,units.kg/units.s),
        "temperature":(216,298,1500,units.K),
        "pressure":(5e4,80*100000,2000*100000,units.Pa)
    },

    "state_definition":FTPx,
    "state_bounds":{
        "flow_mol":(-20000,1,20000,units.mol/units.s),
        "temperature":(216,298,1500,units.K),
        "pressure":(5e4,80*100000,2000*100000,units.Pa)
    },
"""

