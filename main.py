"""This script serves as a template for a Python program that includes a main function."""
import sys
from pprint import pprint

# COMPOUNDS using variable names instead of strings to avoid future typo errors.
# These constants will never change so 'global' is not needed.
SER = 'SER'
THF = 'THF'
NADPP = 'NADP+'
ATP = 'ATP'
CH2THF = 'CH2THF'
GLY = 'GLY'
FORMATE = 'Formate'
ADP = 'ADP'
CHOTHF = 'CHOTHF'
CHPTHF = 'CH+THF'
NADP1P = 'NADP1+'

def mito_ser_thf_to_ch2thf_gly(pool, step, step_timespan):
    """
    This function represents the first step in a biochemical pathway.
    It converts SER and THF into CH2THF and GLY.

    Step 1: [SER] + [THF] <--> [CH2THF] + [GLY]
    """
    ser = pool[SER]  # Retrieve SER from the pool
    thf = pool[THF]  # Retrieve THF from the pool
    ch2thf = pool[CH2THF]  # Retrieve CH2THF from the pool
    gly = pool[GLY]  # Retrieve GLY from the pool
    print(' ', sys._getframe().f_code.co_name)
    print('    ser:', ser, ' thf:', thf, ' ch2thf:', ch2thf, ' gly:', gly)

    # DO: XHESIKA MATH
    ser = ser - 0.1  # Example decrement for SER
    thf = thf - 0.1  # Example decrement for THF
    ch2thf = ch2thf + 0.1  # Example increment for CH2THF
    gly = gly + 0.1  # Example increment for GLY

    pool[SER] = ser
    pool[THF] = thf
    pool[CH2THF] = ch2thf
    pool[GLY] = gly

    print('    ser:', ser, ' thf:', thf, ' ch2thf:', ch2thf, ' gly:', gly)
    return


def mito_ch2thf_nadpp_to_chpthf_nadp(pool, step, step_timespan):
    """
    This function represents the second step in a biochemical pathway.
    It converts CH2THF and NADP+ into CH+THF and NADP1+.

    Step 2a: [CH2THF] + [NADP+] <--> [CH+THF] + [NADP1+]
    """
    ch2thf = pool[CH2THF]  # Retrieve CH2THF from the pool
    nadpp = pool[NADPP]  # Retrieve NADP+ from the pool
    chpthf = pool[CHPTHF]  # Retrieve CH+THF from the pool
    nadp1p = pool[NADP1P]  # Retrieve NADP1+ from the pool
    print(' ', sys._getframe().f_code.co_name)
    print('    ch2thf:', ch2thf, ' nadpp:', nadpp, ' chpthf:', chpthf, ' nadp1p:', nadp1p)

    # DO: XHESIKA MATH
    ch2thf -= 0.1  # Example decrement for CH2THF
    nadpp -= 0.1  # Example decrement for NADP+
    chpthf += 0.1  # Example increment for CH+THF
    nadp1p += 0.1  # Example increment for NADP1+

    pool[CH2THF] = ch2thf
    pool[NADPP] = nadpp
    pool[CHPTHF] = chpthf
    pool[NADP1P] = nadp1p
    print('    ch2thf:', ch2thf, ' nadpp:', nadpp, ' chpthf:', chpthf, ' nadp1p:', nadp1p)
    return


def mito_chpthf_to_chothf(pool, step, step_timespan):
    """
    This function represents the second part of the second step in a biochemical pathway.
    It converts CH+THF into CHOTHF.

    Step 2b: [CH+THF] <--> [CHOTHF]
    """
    chpthf = pool[CHPTHF]  # Retrieve CH+THF from the pool
    chothf = pool[CHOTHF]  # Retrieve CHOTHF from the pool
    print(' ', sys._getframe().f_code.co_name)
    print('    chpthf:', chpthf, ' chothf:', chothf)

    # DO: XHESIKA MATH
    chpthf -= 0.1  # Example decrement for CH+THF
    chothf += 0.1  # Example increment for CHOTHF

    pool[CHPTHF] = chpthf
    pool[CHOTHF] = chothf
    print('    chpthf:', chpthf, ' chothf:', chothf)
    return


def mito_chothf_atp_to_formate_adp(pool, step, step_timespan):
    """
    This function represents the third step in a biochemical pathway.
    It converts CHOTHF and ATP into Formate and ADP.

    Step 3: [CHOTHF] + [ATP] <--> [Formate] + [ADP]
    """
    chothf = pool[CHOTHF]  # Retrieve CHOTHF from the pool
    atp = pool[ATP]  # Retrieve ATP from the pool
    formate = pool[FORMATE]  # Retrieve Formate from the pool
    adp = pool[ADP]  # Retrieve ADP from the pool
    print(' ', sys._getframe().f_code.co_name)
    print('    chothf:', chothf, ' atp:', atp, ' formate:', formate, ' adp:', adp)

    # DO: XHESIKA MATH
    chothf -= 0.1  # Example decrement for CHOTHF
    atp -= 0.1  # Example decrement for ATP
    formate += 0.1  # Example increment for Formate
    adp += 0.1  # Example increment for ADP

    pool[CHOTHF] = chothf
    pool[ATP] = atp
    pool[FORMATE] = formate
    pool[ADP] = adp
    print('    chothf:', chothf, ' atp:', atp, ' formate:', formate, ' adp:', adp)
    return


def main():
    """
    This is the main function that will be executed when the script is run.
    """
    experiment_runtime = 0.1 # Experiment time in seconds
    experiment_timesteps = 4  # Resolution for the simulation
    step_timespan = experiment_runtime / experiment_timesteps

    # Initialize a pool of compounds with their concentrations
    mito_pool = {
        SER : 0.1,
        THF : 0.1,
        NADPP : 0.1,
        ATP : 0.1,
        CH2THF : 0.0,
        GLY : 0.0,
        FORMATE : 0.0,
        ADP : 0.0,
        CHOTHF : 0.0,
        CHPTHF : 0.0,
        NADP1P : 0.0
    }

    print('Initial concentrations:')
    pprint(mito_pool)
    print()

    for step in range(0, experiment_timesteps):
        print('Step:', step, ' Timespan:', step_timespan)
        # Step 1: [SER] + [THF] <--> [CH2THF] + [GLY]
        mito_ser_thf_to_ch2thf_gly(mito_pool, step, step_timespan)

        # Step 2a: [CH2THF] + [NADP+] <--> [CH+THF] + [NADP1+]
        mito_ch2thf_nadpp_to_chpthf_nadp(mito_pool, step, step_timespan)

        # Step 2b: [CH+THF] <--> [CHOTHF]
        mito_chpthf_to_chothf(mito_pool, step, step_timespan)

        # Step 3: [CHOTHF] + [ATP] <--> [Formate] + [ADP]
        mito_chothf_atp_to_formate_adp(mito_pool, step, step_timespan)
        print()

    print('Final concentrations after the experiment:')
    pprint(mito_pool)

# Execute the main function
main()