"""
GlycineCleavageSystem.py

Simplified Glycine Cleavage System model designed to fit the mitochondrial
1-carbon metabolism model structure.

GCS overall reaction:
    GLY + THF + NAD+ --> CH2THF + CO2 + NH3 + NADH + H+

Sequential enzyme steps:
    P-protein: GLY + H_OX --> H_AM + CO2
    T-protein: H_AM + THF --> H_RED + CH2THF + NH3
    L-protein: H_RED + NAD+ --> H_OX + NADH + H+

Units:
    Concentrations: µM
    Rates: µM/h
"""

from ModelParameters import *


class GlycineCleavageSystem:
    """
    Glycine Cleavage System module.

    Calculates instantaneous GCS reaction rates from the supplied
    mitochondrial state.

    This class does not modify metabolite concentrations.

    Stoichiometric concentration changes are handled centrally by
    SerToFormateMmModel.calculate_derivatives().
    """

    def __init__(self):
        pass

    def two_substrate_mm(
            self,
            substrate_1,
            substrate_2,
            vmax,
            km_1,
            km_2
    ):
        """
        Generic two-substrate Michaelis-Menten equation.

        Currently used for GCS P-protein:

            GLY + H_OX --> H_AM + CO2

        Rate:
            v = Vmax
                * S1 / (Km1 + S1)
                * S2 / (Km2 + S2)

        Returns:
            instantaneous rate in µM/h
        """

        if substrate_1 <= 0 or substrate_2 <= 0:
            return 0.0

        return (
                vmax
                * (substrate_1 / (km_1 + substrate_1))
                * (substrate_2 / (km_2 + substrate_2))
        )

    def ordered_bi_bi(self, substrate_a, substrate_b, vmax, km_a, km_b):
        """
        Ordered Bi-Bi equation for T-protein.

        Reaction:
            H_AM + THF --> H_RED + CH2THF + NH3

        A = H_AM
        B = THF

        Approximation:
            KiA ≈ KmA

        Rate:
            v = Vmax * A * B /
                (KiA*KmB + KmB*A + KmA*B + A*B)

        Returns:
            rate in µM/h
        """

        if substrate_a <= 0 or substrate_b <= 0:
            return 0.0

        ki_a = km_a

        denominator = (
            (ki_a * km_b)
            + (km_b * substrate_a)
            + (km_a * substrate_b)
            + (substrate_a * substrate_b)
        )

        return vmax * substrate_a * substrate_b / denominator

    def ping_pong_bi_bi(self, substrate_a, substrate_b, vmax, km_a, km_b):
        """
        Ping-Pong Bi-Bi equation for L-protein.

        Reaction:
            H_RED + NAD+ --> H_OX + NADH + H+

        A = H_RED
        B = NAD+

        Rate:
            v = Vmax * A * B /
                (KmA*B + KmB*A + A*B)

        Returns:
            rate in µM/h
        """

        if substrate_a <= 0 or substrate_b <= 0:
            return 0.0

        denominator = (
            (km_a * substrate_b)
            + (km_b * substrate_a)
            + (substrate_a * substrate_b)
        )

        return vmax * substrate_a * substrate_b / denominator

    def calculate_rates(self, state):
        """
        Calculate all Glycine Cleavage System reaction rates from the
        SAME instantaneous mitochondrial state.

        GCS overall reaction:

            GLY + THF + NAD+
                -->
            CH2THF + CO2 + NH3 + NADH + H+

        Sequential biochemical enzyme reactions:

            P-protein:
                GLY + H_OX --> H_AM + CO2

            T-protein:
                H_AM + THF --> H_RED + CH2THF + NH3

            L-protein:
                H_RED + NAD+ --> H_OX + NADH + H+

        IMPORTANT NUMERICAL PRINCIPLE:

        P-, T-, and L-protein rates are all calculated from the same
        instantaneous mitochondrial state.

        This function does NOT:
            - calculate timestep deltas
            - modify metabolite concentrations
            - consume substrates
            - produce products

        Stoichiometric changes are handled centrally in
        SerToFormateMmModel.calculate_derivatives().

        Returns:
            dictionary containing instantaneous rates in µM/h
        """

        # ==============================================================
        # GCS P-protein
        #
        # GLY + H_OX --> H_AM + CO2
        # ==============================================================

        gly = state[GLY]
        h_ox = state[H_OX]

        rate_p = self.two_substrate_mm(
            substrate_1=gly,
            substrate_2=h_ox,
            vmax=GCS_P_Vmax,
            km_1=GCS_P_KM_GLY,
            km_2=GCS_P_KM_HOX
        )

        # ==============================================================
        # GCS T-protein
        #
        # H_AM + THF --> H_RED + CH2THF + NH3
        # ==============================================================

        h_am = state[H_AM]
        thf = state[THF]

        rate_t = self.ordered_bi_bi(
            substrate_a=h_am,
            substrate_b=thf,
            vmax=GCS_T_Vmax,
            km_a=GCS_T_KM_HAM,
            km_b=GCS_T_KM_THF
        )

        # ==============================================================
        # GCS L-protein
        #
        # H_RED + NAD+ --> H_OX + NADH + H+
        # ==============================================================

        h_red = state[H_RED]
        nad_plus = state[NAD_PLUS]

        rate_l = self.ping_pong_bi_bi(
            substrate_a=h_red,
            substrate_b=nad_plus,
            vmax=GCS_L_Vmax,
            km_a=GCS_L_KM_HRED,
            km_b=GCS_L_KM_NAD
        )

        # ==============================================================
        # Return rates only
        # ==============================================================

        return {
            'gcs_p_rate': rate_p,
            'gcs_t_rate': rate_t,
            'gcs_l_rate': rate_l
        }


