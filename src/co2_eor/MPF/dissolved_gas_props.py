import pyomo.environ as pyo
from idaes.core.util.misc import set_param_from_config
from idaes.models.properties.modular_properties.pure import NIST


# -----------------------------------------------------------------------------
# custom enthalpy calc
class custom_calc(object):

    class enth_mol_liq_comp(object):
        @staticmethod
        def build_parameters(cobj):
            #if not hasattr(cobj, "cp_mol_liq_comp_coeff"):
            #    custom_calc.cp_mol_liq_comp.build_parameters(cobj)
            units = cobj.parent_block().get_metadata().derived_units

            #get enthalpy of dissolution into phase
            cobj.enth_diss = pyo.Var(doc="heat of dissolution",units=units.ENERGY_MOLE)
            set_param_from_config(cobj, param="enth_diss")

            if cobj.parent_block().config.include_enthalpy_of_formation:
                #units = cobj.parent_block().get_metadata().derived_units

                cobj.enth_mol_form_liq_comp_ref = pyo.Var(
                    doc="Liquid phase molar heat of formation @ Tref",
                    units=units.ENERGY_MOLE,
                )
                set_param_from_config(cobj, param="enth_mol_form_liq_comp_ref")

        @staticmethod
        def return_expression(b, cobj, T):
            # Specific enthalpy
            units = b.params.get_metadata().derived_units
            Tr = b.params.temperature_ref



            h_form = (
                cobj.enth_mol_form_liq_comp_ref
                if b.params.config.include_enthalpy_of_formation
                else 0 * units.ENERGY_MOLE
            )

            return NIST.enth_mol_ig_comp.return_expression(b, cobj,T) + cobj.enth_diss + h_form

    class dens_mol_liq_comp(object):
        @staticmethod
        def build_parameters(cobj):
            units = cobj.parent_block().get_metadata().derived_units
            cobj.dens_mol_liq_comp_coeff = pyo.Var(
                doc="Parameter for liquid phase molar density",
                units=units.DENSITY_MOLE,
            )
            set_param_from_config(cobj, param="dens_mol_liq_comp_coeff")

        @staticmethod
        def return_expression(b, cobj, T):
            # Molar density
            return cobj.dens_mol_liq_comp_coeff