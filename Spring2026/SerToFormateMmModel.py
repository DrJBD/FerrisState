"""
Report which is a ReportingTable should be in the globals when running this.
"""
from pprint import pprint
from ModelParameters import *
from GlycineCleavageSystem import GlycineCleavageSystem
from ReportingTable import report
import numpy as np
from scipy.integrate import solve_ivp


class SerToFormateMmModel:

    def __init__(self, mito_pool, experiment_runtime_in_sec, experiment_timesteps):
        self.mito_pool = mito_pool
        self.experiment_runtime_in_sec = experiment_runtime_in_sec
        self.experiment_timesteps = experiment_timesteps
        self.gcs = GlycineCleavageSystem()
        self.ran = False  # Flag to indicate if the model has been run

    def get_state_compounds(self):
        """
        Fixed metabolite ordering used by the ODE solver.

        This ordering must remain identical when converting
        dictionary states to solver vectors and back.
        """

        return STATE_COMPOUNDS

    def state_to_vector(self, state):
        """
        Convert metabolite dictionary to numerical solver vector.
        """

        compounds = self.get_state_compounds()

        return [
            state[compound]
            for compound in compounds
        ]

    def vector_to_state(self, y):
        """
        Convert numerical solver vector back to metabolite dictionary.
        """

        compounds = self.get_state_compounds()

        return {
            compound: concentration
            for compound, concentration in zip(compounds, y)
        }

    def ode_system(self, _t, y):
        """
        ODE system passed to scipy.integrate.solve_ivp.

        For a given instantaneous state y:

            1. Convert vector -> metabolite dictionary.
            2. Calculate every reaction rate from that same state.
            3. Calculate total metabolite derivatives.
            4. Return derivatives as a numerical vector.

        Time t is supplied by solve_ivp.
        Reaction rates and derivatives are in µM/h,
        so t is expressed in hours.
        """

        state = self.vector_to_state(y)

        rates = self.calculate_all_rates(state)

        derivatives = self.calculate_derivatives(
            state,
            rates
        )

        compounds = self.get_state_compounds()

        return [
            derivatives[compound]
            for compound in compounds
        ]


    def km_ser_app_from_thf(self, thf):
        """
        Apparent Km for serine as a function of THF concentration.
        Based on Fig. 2F.

        thf: [THF] in µM
        returns: Km_ser_app in µM
        """
        return (
                Step1_ALPHA_KM_SER
                + (Step1_KM_SER_MAX - Step1_ALPHA_KM_SER) * Step1_K05_THF / (Step1_K05_THF + thf)
        )

    def shmt2_forward_inhibited(self, thf, ser, Step1_Vmaxf, Step1_aKm_thf, Step1_Ki_thf, Step1_Km_ser_app, Step1_Ki_ser):
        """
        Calculates the rate of the SHMT2 forward reaction with substrate inhibition.

        SER + THF --> CH2THF + GLY
        """
        thf_term = thf / (Step1_aKm_thf + thf * (1 + thf / Step1_Ki_thf))
        ser_term = ser / (Step1_Km_ser_app + ser * (1 + ser / Step1_Ki_ser))
        return Step1_Vmaxf * thf_term * ser_term

    def shmt2_reverse_rate(self, gly, ch2thf, Step1_Vmaxr, Step1_aKm_gly, Step1_Km_ch2thf):
        """
        Calculates the rate of the SHMT2 reverse reaction with Michaelis-Menten kinetics.

        GLY + CH2THF --> SER + THF
        """
        gly_term = gly / (Step1_aKm_gly + gly)
        ch2thf_term = ch2thf / (Step1_Km_ch2thf + ch2thf)
        return Step1_Vmaxr * gly_term * ch2thf_term

    def calculate_shmt2_rates(self, state):
        """
        Calculate SHMT2 forward and reverse reaction rates from the same
        instantaneous mitochondrial state.

        Reaction:

            SER + THF <--> CH2THF + GLY

        Forward:
            SER + THF --> CH2THF + GLY

        Reverse:
            GLY + CH2THF --> SER + THF

        IMPORTANT:
        This function calculates rates only.

        It does NOT:
            - calculate delta
            - modify the metabolite pool
            - consume substrates
            - produce products

        Pool changes will later be handled centrally by
        calculate_derivatives() and the numerical integrator.

        Returns:
            dictionary containing:
                SHMT2 forward rate
                SHMT2 reverse rate
                SHMT2 net rate

        Units:
            µM/h
        """

        # --------------------------------------------------------------
        # Read concentrations from ONE instantaneous state
        # --------------------------------------------------------------

        ser = state[SER]
        thf = state[THF]
        gly = state[GLY]
        ch2thf = state[CH2THF]

        # --------------------------------------------------------------
        # Calculate apparent Km for serine from instantaneous THF
        # --------------------------------------------------------------

        km_ser_app = self.km_ser_app_from_thf(thf)

        # --------------------------------------------------------------
        # SHMT2 forward rate
        #
        # SER + THF --> CH2THF + GLY
        # --------------------------------------------------------------

        rate_forward = self.shmt2_forward_inhibited(
            thf=thf,
            ser=ser,
            Step1_Vmaxf=Step1_Vmax_f,
            Step1_aKm_thf=Step1_aKm_thf,
            Step1_Ki_thf=Step1_Ki_thf,
            Step1_Km_ser_app=km_ser_app,
            Step1_Ki_ser=Step1_Ki_ser
        )

        # --------------------------------------------------------------
        # SHMT2 reverse rate
        #
        # GLY + CH2THF --> SER + THF
        # --------------------------------------------------------------

        rate_reverse = self.shmt2_reverse_rate(
            gly=gly,
            ch2thf=ch2thf,
            Step1_Vmaxr=Step1_Vmax_r,
            Step1_aKm_gly=Step1_aKm_gly,
            Step1_Km_ch2thf=Step1_Km_ch2thf
        )

        # --------------------------------------------------------------
        # Net SHMT2 rate
        #
        # Positive = net forward
        # Negative = net reverse
        # --------------------------------------------------------------

        net_rate = (
                rate_forward
                - rate_reverse
        )

        return {
            'shmt2_forward_rate': rate_forward,
            'shmt2_reverse_rate': rate_reverse,
            'shmt2_net_rate': net_rate
        }


    # ==================================================================
    # STEP 2a
    #
    # MTHFD2 and MTHFD2L dehydrogenase activities
    # ==================================================================
    #
    # Four distinct reaction fluxes are modeled:
    #
    # 1. MTHFD2 + NAD+
    #
    #       CH2THF + NAD+ --> CH+THF + NADH
    #
    # 2. MTHFD2 + NADP+
    #
    #       CH2THF + NADP+ --> CH+THF + NADPH
    #
    # 3. MTHFD2L + NAD+
    #
    #       CH2THF + NAD+ --> CH+THF + NADH
    #
    # 4. MTHFD2L + NADP+
    #
    #       CH2THF + NADP+ --> CH+THF + NADPH
    #
    #
    # MTHFD2 and MTHFD2L are independent enzyme populations.
    #
    # However, NAD+ and NADP+ are treated as alternative competing
    # cofactors within each enzyme.
    #
    # All four rates are calculated from the SAME instantaneous
    # mitochondrial state.
    #
    # No metabolite concentrations are modified inside these functions.
    # ==================================================================

    def competitive_cofactor_rate(
            self,
            ch2thf,
            cofactor,
            competing_cofactor,
            vmax,
            km_ch2thf,
            km_cofactor,
            ki_competing_cofactor
    ):

        """
        Calculate one cofactor-specific MTHFD2/MTHFD2L reaction rate.

        General form:

                               [CH2THF]
            v = Vmax * ---------------------------
                         Km_CH2THF + [CH2THF]

                               [cofactor]
                     * ----------------------------------------
                       Km_cofactor *
                       (1 + [competitor] / Ki_competitor)
                       + [cofactor]


        Example for an NAD+-dependent branch:

                               [CH2THF]
            v_NAD = Vmax * -------------------------
                           Km_CH2THF + [CH2THF]

                                  [NAD+]
                        * ------------------------------------
                          Km_NAD *
                          (1 + [NADP+] / Ki_NADP)
                          + [NAD+]


        Biological/modeling interpretation
        -----------------------------------

        NAD+ and NADP+ are treated as alternative competing cofactors
        for the same enzyme population.

        MTHFD2 and MTHFD2L themselves are NOT treated as competitive
        inhibitors of one another. They are separate enzymes.

        IMPORTANT:
        Experimental mixed-cofactor Ki values are currently unavailable.

        Therefore, the baseline model uses:

            Ki_competing_cofactor ≈ Km_competing_cofactor

        Km and Ki remain separate parameters in the code because they
        are not fundamentally equivalent.

        Units:
            concentrations: µM
            Vmax: µM/h
            returned rate: µM/h
        """

        if ch2thf <= 0 or cofactor <= 0:
            return 0.0

        # Prevent a small numerical negative concentration from entering
        # the kinetic equation.
        ch2thf = max(ch2thf, 0.0)
        cofactor = max(cofactor, 0.0)
        competing_cofactor = max(competing_cofactor, 0.0)

        # --------------------------------------------------------------
        # CH2THF saturation
        # --------------------------------------------------------------

        ch2thf_term = (
                ch2thf
                / (
                        km_ch2thf
                        + ch2thf
                )
        )

        # --------------------------------------------------------------
        # Cofactor term including competition from the alternative
        # redox cofactor
        # --------------------------------------------------------------

        cofactor_term = (
                cofactor
                / (
                        km_cofactor
                        * (
                                1.0
                                + competing_cofactor
                                / ki_competing_cofactor
                        )
                        + cofactor
                )
        )

        # --------------------------------------------------------------
        # Final instantaneous rate
        # --------------------------------------------------------------

        rate = (
                vmax
                * ch2thf_term
                * cofactor_term
        )

        return rate

    def calculate_mthfd2_mthfd2l_rates(self, state):
        """
        Calculate all four MTHFD2/MTHFD2L Step 2a fluxes.

        IMPORTANT:
        Every rate is calculated from the SAME instantaneous state.

        The state is not modified anywhere in this function.

        Returns:
            dictionary containing four separate instantaneous rates,
            all in µM/h.
        """

        # --------------------------------------------------------------
        # Read ONE instantaneous mitochondrial state
        # --------------------------------------------------------------

        ch2thf = state[CH2THF]

        nad = state[NAD_PLUS]
        nadp = state[NADPP]

        # ==============================================================
        # 1. MTHFD2 + NAD+
        #
        # CH2THF + NAD+ --> CH+THF + NADH
        # ==============================================================

        mthfd2_nad_rate = self.competitive_cofactor_rate(

            ch2thf=ch2thf,

            cofactor=nad,
            competing_cofactor=nadp,

            vmax=MTHFD2_NAD_VMAX,

            km_ch2thf=MTHFD2_NAD_KM_CH2THF,

            km_cofactor=MTHFD2_NAD_KM_NAD,

            ki_competing_cofactor=MTHFD2_NAD_KI_NADP
        )

        # ==============================================================
        # 2. MTHFD2 + NADP+
        #
        # CH2THF + NADP+ --> CH+THF + NADPH
        # ==============================================================

        mthfd2_nadp_rate = self.competitive_cofactor_rate(

            ch2thf=ch2thf,

            cofactor=nadp,
            competing_cofactor=nad,

            vmax=MTHFD2_NADP_VMAX,

            km_ch2thf=MTHFD2_NADP_KM_CH2THF,

            km_cofactor=MTHFD2_NADP_KM_NADP,

            ki_competing_cofactor=MTHFD2_NADP_KI_NAD
        )

        # ==============================================================
        # 3. MTHFD2L + NAD+
        #
        # CH2THF + NAD+ --> CH+THF + NADH
        # ==============================================================

        mthfd2l_nad_rate = self.competitive_cofactor_rate(

            ch2thf=ch2thf,

            cofactor=nad,
            competing_cofactor=nadp,

            vmax=MTHFD2L_NAD_VMAX,

            km_ch2thf=MTHFD2L_NAD_KM_CH2THF,

            km_cofactor=MTHFD2L_NAD_KM_NAD,

            ki_competing_cofactor=MTHFD2L_NAD_KI_NADP
        )

        # ==============================================================
        # 4. MTHFD2L + NADP+
        #
        # CH2THF + NADP+ --> CH+THF + NADPH
        # ==============================================================

        mthfd2l_nadp_rate = self.competitive_cofactor_rate(

            ch2thf=ch2thf,

            cofactor=nadp,
            competing_cofactor=nad,

            vmax=MTHFD2L_NADP_VMAX,

            km_ch2thf=MTHFD2L_NADP_KM_CH2THF,

            km_cofactor=MTHFD2L_NADP_KM_NADP,

            ki_competing_cofactor=MTHFD2L_NADP_KI_NAD
        )

        # --------------------------------------------------------------
        # Return all four fluxes separately.
        #
        # Do NOT collapse them into one rate here because we want
        # individual reporting of:
        #
        #   - enzyme identity
        #   - cofactor identity
        #   - NADH-producing flux
        #   - NADPH-producing flux
        # --------------------------------------------------------------

        return {

            'mthfd2_nad_rate':
                mthfd2_nad_rate,

            'mthfd2_nadp_rate':
                mthfd2_nadp_rate,

            'mthfd2l_nad_rate':
                mthfd2l_nad_rate,

            'mthfd2l_nadp_rate':
                mthfd2l_nadp_rate
        }

    def irreversible_two_substrate_mm(self, substrate_1, substrate_2, vmax, km_1, km_2):
        """
        Calculates an irreversible two-substrate Michaelis-Menten rate.

        General reaction:
            S1 + S2 --> products

        Rate:
            v = Vmax * [S1]/(Km1 + [S1]) * [S2]/(Km2 + [S2])

        Units:
            substrate concentrations: µM
            Vmax: µM/h
            returned rate: µM/h
        """

        if substrate_1 <= 0 or substrate_2 <= 0:
            return 0.0

        return (
                vmax
                * (substrate_1 / (km_1 + substrate_1))
                * (substrate_2 / (km_2 + substrate_2))
        )

    def irreversible_single_substrate_mm(self, substrate, vmax, km):
        """
        Calculates an irreversible single-substrate Michaelis-Menten rate.

        General reaction:
            S --> products

        Rate:
            v = Vmax * [S]/(Km + [S])

        Used for Step 2b because H2O is treated as constant and non-limiting.
        """

        if substrate <= 0:
            return 0.0

        return vmax * (substrate / (km + substrate))

    def calculate_all_rates(self, state):
        """
        Calculate every pathway reaction rate from the same
        instantaneous mitochondrial state.

        No concentrations are modified here.
        """

        rates = {}

        # Step 1 - SHMT2
        rates.update(
            self.calculate_shmt2_rates(state)
        )

        # GCS
        rates.update(
            self.gcs.calculate_rates(state)
        )

        # Step 2a - MTHFD2 / MTHFD2L
        rates.update(
            self.calculate_mthfd2_mthfd2l_rates(state)
        )

        # Step 2b
        rates['step2b_rate'] = self.calculate_step2b_rate(state)

        # Step 3
        rates['step3_rate'] = self.calculate_step3_rate(state)

        # Step 4
        rates['step4_rate'] = self.calculate_step4_rate(state)

        return rates

    def calculate_derivatives(self, state, rates):
        """
        Convert instantaneous reaction rates into concentration derivatives.

        All derivatives are calculated from the same instantaneous state.

        No metabolite pool is updated inside this function.

        Returns:
            Dictionary containing d[metabolite]/dt values in µM/h.
        """

        # Start every metabolite derivative at zero.
        derivatives = {
            compound: 0.0
            for compound in state
        }

        # ==================================================================
        # Step 1
        # SHMT2
        #
        # SER + THF <--> CH2THF + GLY
        # ==================================================================

        shmt2_forward = rates['shmt2_forward_rate']
        shmt2_reverse = rates['shmt2_reverse_rate']

        # SER
        derivatives[SER] += (
                -shmt2_forward
                + shmt2_reverse
        )

        # THF
        derivatives[THF] += (
                -shmt2_forward
                + shmt2_reverse
        )

        # CH2THF
        derivatives[CH2THF] += (
                shmt2_forward
                - shmt2_reverse
        )

        # GLY
        derivatives[GLY] += (
                shmt2_forward
                - shmt2_reverse
        )

        # ==================================================================
        # Step 2a
        # MTHFD2 / MTHFD2L dehydrogenase reactions
        # ==================================================================

        mthfd2_nad = rates['mthfd2_nad_rate']
        mthfd2_nadp = rates['mthfd2_nadp_rate']

        mthfd2l_nad = rates['mthfd2l_nad_rate']
        mthfd2l_nadp = rates['mthfd2l_nadp_rate']

        # --------------------------------------------------------------
        # Total NAD+-dependent flux
        # --------------------------------------------------------------

        total_nad_flux = (
                mthfd2_nad
                + mthfd2l_nad
        )

        # --------------------------------------------------------------
        # Total NADP+-dependent flux
        # --------------------------------------------------------------

        total_nadp_flux = (
                mthfd2_nadp
                + mthfd2l_nadp
        )

        # --------------------------------------------------------------
        # Total Step 2a flux
        # --------------------------------------------------------------

        total_step2a_flux = (
                total_nad_flux
                + total_nadp_flux
        )

        # CH2THF consumed by all four reactions
        derivatives[CH2THF] -= total_step2a_flux

        # CH+THF produced by all four reactions
        derivatives[CHPTHF] += total_step2a_flux

        # NAD+-dependent branches
        derivatives[NAD_PLUS] -= total_nad_flux
        derivatives[NADH] += total_nad_flux

        # NADP+-dependent branches
        derivatives[NADPP] -= total_nadp_flux
        derivatives[NADPH] += total_nadp_flux

        # ==================================================================
        # Glycine Cleavage System
        # ==================================================================

        # ==================================================================
        # GCS P-protein
        #
        # GLY + H_OX --> H_AM + CO2
        # ==================================================================

        gcs_p = rates['gcs_p_rate']

        derivatives[GLY] -= gcs_p
        derivatives[H_OX] -= gcs_p

        derivatives[H_AM] += gcs_p
        derivatives[CO2] += gcs_p

        # ==================================================================
        # GCS T-protein
        #
        # H_AM + THF --> H_RED + CH2THF + NH3
        # ==================================================================

        gcs_t = rates['gcs_t_rate']

        derivatives[H_AM] -= gcs_t
        derivatives[THF] -= gcs_t

        derivatives[H_RED] += gcs_t
        derivatives[CH2THF] += gcs_t
        derivatives[NH3] += gcs_t

        # ==================================================================
        # GCS L-protein
        #
        # H_RED + NAD+ --> H_OX + NADH + H+
        # ==================================================================

        gcs_l = rates['gcs_l_rate']

        derivatives[H_RED] -= gcs_l
        derivatives[NAD_PLUS] -= gcs_l

        derivatives[H_OX] += gcs_l
        derivatives[NADH] += gcs_l

        # H+ handling is temporary.
        # We will review whether H+ should remain a dynamic pool.
        derivatives[H_PLUS] += gcs_l

        # =================================================================
        # Step 2b
        # CH+THF --> CHOTHF + H+
        # ==================================================================

        step2b = rates['step2b_rate']

        derivatives[CHPTHF] -= step2b
        derivatives[CHOTHF] += step2b

        # Temporary proton bookkeeping
        derivatives[H_PLUS] += step2b

        # ==================================================================
        # Step 3
        # CHOTHF + ADP --> Formate + THF + ATP
        # ==================================================================

        step3 = rates['step3_rate']

        derivatives[CHOTHF] -= step3
        derivatives[ADP] -= step3

        derivatives[FORMATE] += step3
        derivatives[THF] += step3
        derivatives[ATP] += step3

        # ==================================================================
        # Step 4
        # CHOTHF --> Formate + THF
        # ==================================================================

        step4 = rates['step4_rate']

        derivatives[CHOTHF] -= step4
        derivatives[FORMATE] += step4
        derivatives[THF] += step4

        return derivatives

    def calculate_step2b_rate(self, state):
        """
        Step 2b:
            CH+THF --> CHOTHF + H+

        MTHFD2/MTHFD2L cyclohydrolase activity.

        Water is treated as constant and non-limiting, so the rate law
        depends only on CH+THF.

        Returns:
            instantaneous rate in µM/h
        """

        chpthf = state[CHPTHF]

        rate = self.irreversible_single_substrate_mm(
            substrate=chpthf,
            vmax=Step2b_Vmax_f,
            km=Step2b_Km_chpthf
        )

        return rate

    def calculate_step3_rate(self, state):
        """
        Step 3:
            CHOTHF + ADP --> Formate + THF + ATP

        Simplified MTHFD1L formate-producing direction.

        Pi is not currently represented explicitly.

        Returns:
            instantaneous rate in µM/h
        """

        chothf = state[CHOTHF]
        adp = state[ADP]

        rate = self.irreversible_two_substrate_mm(
            substrate_1=chothf,
            substrate_2=adp,
            vmax=Step3_Vmax_f,
            km_1=Step3_Km_chothf,
            km_2=Step3_Km_adp
        )

        return rate

    def calculate_step4_rate(self, state):
        """
        Step 4:
            CHOTHF + H2O --> Formate + THF

        ALDH1L2/FDH hydrolase activity.

        Water is treated as constant and non-limiting.

        Returns:
            instantaneous rate in µM/h
        """

        chothf = state[CHOTHF]

        rate = self.irreversible_single_substrate_mm(
            substrate=chothf,
            vmax=Step4_Vmax,
            km=Step4_Km_chothf
        )

        return rate

    def ser_to_formate_mm_equation(
            self,
            pool,
            experiment_runtime_in_sec,
            experiment_timesteps
    ):
        """
        Integrate the mitochondrial serine-to-formate pathway using
        scipy.integrate.solve_ivp.

        All reaction rates are evaluated from the instantaneous state
        supplied by the ODE solver.

        Internal solver timesteps are adaptive and are independent of
        reporting times.
        """

        print('Initial concentrations:')
        pprint(pool)
        print()

        # --------------------------------------------------------------
        # Convert experiment time to hours
        #
        # All kinetic rates are expressed in µM/h.
        # --------------------------------------------------------------

        experiment_runtime_hr = (
                experiment_runtime_in_sec
                / 3600.0
        )

        # --------------------------------------------------------------
        # Initial solver state
        # --------------------------------------------------------------

        y0 = self.state_to_vector(pool)

        # --------------------------------------------------------------
        # Reporting times
        #
        # These are NOT solver integration steps.
        #
        # solve_ivp chooses its own internal adaptive timesteps.
        # t_eval only specifies where we want results returned.
        # --------------------------------------------------------------

        t_eval = np.linspace(
            0.0,
            experiment_runtime_hr,
            experiment_timesteps + 1
        )

        # --------------------------------------------------------------
        # Solve ODE system
        # --------------------------------------------------------------

        solution = solve_ivp(
            fun=self.ode_system,
            t_span=(
                0.0,
                experiment_runtime_hr
            ),
            y0=y0,
            method='RK45',
            t_eval=t_eval,
            rtol=1e-6,
            atol=1e-9
        )

        # --------------------------------------------------------------
        # Check solver success
        # --------------------------------------------------------------

        if not solution.success:
            raise RuntimeError(
                f'ODE solver failed: {solution.message}'
            )

        compounds = self.get_state_compounds()

        # --------------------------------------------------------------
        # Report each requested time point
        # --------------------------------------------------------------

        for index, time_hr in enumerate(solution.t):
            state = {
                compound: solution.y[i, index]
                for i, compound in enumerate(compounds)
            }

            rates = self.calculate_all_rates(state)

            report_data = {
                **state,
                **rates,
                'time_hr': time_hr,
                'time_sec': time_hr * 3600.0
            }

            report.add_step(
                index,
                report_data
            )

        # --------------------------------------------------------------
        # Final state
        # --------------------------------------------------------------

        final_state = {
            compound: solution.y[i, -1]
            for i, compound in enumerate(compounds)
        }

        self.mito_pool.update(final_state)

        print('Final concentrations after the experiment:')
        pprint(final_state)

        return final_state

    def run(self):
        if not self.ran:
            self.ran = True
            return self.ser_to_formate_mm_equation(self.mito_pool, self.experiment_runtime_in_sec, self.experiment_timesteps)
        else:
            print("Model has already been run. Returning the previous results.")
            return self.mito_pool
