import sys
from pprint import pprint
from ModelParameters import *
from GlycineCleavageSystem import GlycineCleavageSystem

class SerToFormateMmModel:

    def __init__(self, mito_pool, experiment_runtime_in_sec, experiment_timesteps):
        self.mito_pool = mito_pool
        self.experiment_runtime_in_sec = experiment_runtime_in_sec
        self.experiment_timesteps = experiment_timesteps
        self.ran = False  # Flag to indicate if the model has been run

    def two_substrate_mm(self, A, B, Km_A, Km_B, Vmax):
        """
        This function computes the rate for a two-substrate reaction
        A: Concentration of substrate A
        B: Concentration of substrate B
        Vmax: Maximum rate of reaction.
        Km_A, Km_B: Taken from literature. Constants. Dictate the concentration of substrates as half-Vmax is reached.
        """
        return (Vmax * A * B) / (Km_A * B + Km_B * A + A * B)


    def mito_ser_thf_to_ch2thf_gly(self, pool, step, step_timespan):
        """
        This function represents the first step in a biochemical pathway.
        It converts SER and THF into CH2THF and GLY.

        Step 1: [SER] + [THF] <--> [CH2THF] + [GLY]
        """
        print(' ', sys._getframe().f_code.co_name)
        print('    ser:', pool[SER], ' thf:', pool[THF], ' ch2thf:', pool[CH2THF], ' gly:', pool[GLY])

        A = pool[SER]  # input  # Retrieve SER from the pool (uM)
        B = pool[THF]  # Retrieve THF from the pool (uM)
        # pool[Gly]          #Initial Gly concentration (uM)
        # pool[CH2THF]       #Initial 5-10-CH2-THF conc (uM)

        # Constants for enzyme kinetics
        Vmax_1 = 11440  # uM/hour
        Km_ser = 600  # uM
        Km_thf = 50  # uM

        # Calculate rate using two substrate MM Equation
        rate_1 = self.two_substrate_mm(A, B, Km_ser, Km_thf, Vmax_1)
        delta_1 = rate_1 * step_timespan

        pool[SER] -= delta_1
        pool[THF] -= delta_1
        pool[CH2THF] += delta_1
        pool[GLY] += delta_1

        Cleavage = GlycineCleavageSystem()
        Cleavage.mito_gcs_step(pool, step_timespan)

        print('    ser:', pool[SER], ' thf:', pool[THF], ' ch2thf:', pool[CH2THF], ' gly:', pool[GLY])
        return


    def mito_ch2thf_nadpp_to_chpthf_nadp(self, pool, step, step_timespan):
        """
        This function represents the second step in a biochemical pathway.
        It converts CH2THF and NADP+ into CH+THF and NADPH.

        Step 2a: [CH2THF] + [NADP+] <--> [CH+THF] + [NADPH]
        """
        print(' ', sys._getframe().f_code.co_name)
        print('    ch2thf:', pool[CH2THF], ' nadpp:', pool[NADPP], ' chpthf:', pool[CHPTHF], ' nadph:', pool[NADPH])

        A = pool[CH2THF]  # uM (This pool should resemble the updated concentrations (pool[CH2THF] += delta_1)
        B = pool[NADPP]  # uM Approximation
        # pool[CHPTHF]    #Initial intermediate concentration (uM)
        # pool[NADPH]      #uM Approximation

        # Constants for enzyme kinetics
        Vmax_2a = 15000  # unknown for the moment. Catalytic efficiency "kcat/Km" is 0.067 (1/sM). The given value is an approximation
        Km_ch2thf = 50  # uM
        Km_nadpp = 537  # uM

        # Calculate rate using two substrate MM equation
        rate_2a = self.two_substrate_mm(A, B, Vmax_2a, Km_ch2thf, Km_nadpp)
        delta_2a = rate_2a * step_timespan

        # New concentrations as reaction proceeds
        pool[CH2THF] -= delta_2a
        pool[NADPP] -= delta_2a
        pool[CHPTHF] += delta_2a
        pool[NADPH] += delta_2a

        print('    ch2thf:', pool[CH2THF], ' nadpp:', pool[NADPP], ' chpthf:', pool[CHPTHF], ' nadph:', pool[NADPH])
        return


    def mito_chpthf_to_chothf(self, pool, step, step_timespan):
        """
        This function represents the second part of the second step in a biochemical pathway.
        It converts CH+THF into CHOTHF.

        Step 2b: [CH+THF] + [H2O] <--> [CHOTHF] + [H+]
        """
        print(' ', sys._getframe().f_code.co_name)
        print('    chpthf:', pool[CHPTHF], ' h2o:', pool[H2O], ' chothf:', pool[CHOTHF], ' h_plus:', pool[H_PLUS])

        Vmax_2b = 20000  # uM/h
        Km_chpthf = 50  # uM Approximation
        Km_h2o = 0.001  # uM Approximation

        A = pool[CHPTHF]  # updated values from "pool[CHPTHF] += delta_2a"
        B = pool[H2O]  # Excess (uM). Approximation
        # pool[CHOTHF]        #Initial concentration
        # pool[H_PLUS]        #Excess (uM). Approximation.

        # Calculate rate using two substrate MM equation
        rate_2b = self.two_substrate_mm(A, B, Vmax_2b, Km_chpthf, Km_h2o)
        delta_2b = rate_2b * step_timespan

        # New concentrations as reaction proceeds
        pool[CHPTHF] -= delta_2b
        pool[H2O] -= delta_2b
        pool[CHOTHF] += delta_2b
        pool[H_PLUS] += delta_2b

        print('    chpthf:', pool[CHPTHF], ' h2o:', pool[H2O], ' chothf:', pool[CHOTHF], ' h_plus:', pool[H_PLUS])
        return


    def mito_chothf_atp_to_formate_adp(self, pool, step, step_timespan):
        """
        This function represents the third step in a biochemical pathway.
        It converts CHOTHF and ATP into Formate and ADP.

        Step 3: [CHOTHF] + [ATP] <--> [Formate] + [ADP]
        """
        print(' ', sys._getframe().f_code.co_name)
        print('    chothf:', pool[CHOTHF], ' atp:', pool[ATP], ' formate:', pool[FORMATE], ' adp:', pool[ADP])

        Vmax_3 = 11440  # uM/h Approximation.
        Km_chothf = 50  # uM Approximation
        Km_atp = 40  # uM

        A = pool[CHOTHF]  # Updated (uM)
        B = pool[ATP]  # Excess (uM) Approximation
        # pool[FORMATE]          #(uM) If we're trying to isolate the net production of formate
        # pool[ADP]              #Excess (uM) Approximation

        # Calculate rate using two substrate MM equation
        rate_3 = self.two_substrate_mm(A, B, Vmax_3, Km_chothf, Km_atp)
        delta_3 = rate_3 * step_timespan

        # New concentrations as reaction proceeds
        pool[CHOTHF] -= delta_3
        pool[ATP] -= delta_3
        pool[FORMATE] += delta_3  # Outcome
        pool[ADP] += delta_3

        print('    chothf:', pool[CHOTHF], ' atp:', pool[ATP], ' formate:', pool[FORMATE], ' adp:', pool[ADP])
        return


    def ser_to_formate_mm_equation(self, pool, experiment_runtime_in_sec, experiment_timesteps):
        """
        This function calculates the concentration of Formate produced from a given concentration of SER

        Returns the updated pool of concentrations including the Formate produced.
        """
        print('Initial concentrations:')
        pprint(pool)
        print()

        step_timespan = experiment_runtime_in_sec / experiment_timesteps
        for step in range(0, experiment_timesteps):
            print('Step:', step, ' Timespan:', step_timespan)
            # Step 1: [SER] + [THF] <--> [CH2THF] + [GLY]
            self.mito_ser_thf_to_ch2thf_gly(pool, step, step_timespan)

            # Step 2a: [CH2THF] + [NADP+] <--> [CH+THF] + [NADPH]
            self.mito_ch2thf_nadpp_to_chpthf_nadp(pool, step, step_timespan)

            # Step 2b: [CH+THF] <--> [CHOTHF]
            self.mito_chpthf_to_chothf(pool, step, step_timespan)

            # Step 3: [CHOTHF] + [ATP] <--> [Formate] + [ADP]
            self.mito_chothf_atp_to_formate_adp(pool, step, step_timespan)
            print()

        print('Final concentrations after the experiment:')
        pprint(pool)
        return pool  # Return the concentration of Formate produced


    def run(self):
        if not self.ran:
            self.ran = True
            return self.ser_to_formate_mm_equation(self.mito_pool, self.experiment_runtime_in_sec, self.experiment_timesteps)
        else:
            print("Model has already been run. Returning the previous results.")
            return self.mito_pool
