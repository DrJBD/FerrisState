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
H_PLUS = 'H+'
FORMATE = 'Formate'
ADP = 'ADP'
CHOTHF = 'CHOTHF'
CHPTHF = 'CH+THF'
NADPH = 'NADPH'
H2O = 'H2O'

def two_substrate_mm(A, B, Km_A, Km_B, Vmax):
    """
    This function computes the rate for a two-substrate reaction
    A: Concentration of substrate A
    B: Concentration of substrate B
    Vmax: Maximum rate of reaction.
    Km_A, Km_B: Taken from literature. Constants. Dictate the concentration of substrates as half-Vmax is reached.
    """
    return (Vmax * A * B) / (Km_A * B + Km_B * A + A * B)

def mito_ser_thf_to_ch2thf_gly(pool, step, step_timespan):
    """
    This function represents the first step in a biochemical pathway.
    It converts SER and THF into CH2THF and GLY.

    Step 1: [SER] + [THF] <--> [CH2THF] + [GLY]
    """
    print(' ', sys._getframe().f_code.co_name)
    print('    ser:', pool[SER], ' thf:', pool[THF], ' ch2thf:', pool[CH2THF], ' gly:', pool[GLY])

    A = pool[SER]  # Retrieve SER from the pool
    B = pool[THF]  # Retrieve THF from the pool

    #Constants for enzyme kinetics
    Vmax_1 = 1.0
    Km_ser = 0.05
    Km_thf = 0.05

    #Caculate rate using two substrate MM Equation
    rate_1 = two_substrate_mm(A, B, Km_ser, Km_thf, Vmax_1)
    delta_1 = rate_1 * step_timespan

    pool[SER] -= delta_1
    pool[THF] -= delta_1
    pool[CH2THF] += delta_1
    pool[GLY] += delta_1

    print('    ser:', pool[SER], ' thf:', pool[THF], ' ch2thf:', pool[CH2THF], ' gly:', pool[GLY])
    return


def mito_ch2thf_nadpp_to_chpthf_nadp(pool, step, step_timespan):
    """
    This function represents the second step in a biochemical pathway.
    It converts CH2THF and NADP+ into CH+THF and NADPH.

    Step 2a: [CH2THF] + [NADP+] <--> [CH+THF] + [NADPH] 
    """
    print(' ', sys._getframe().f_code.co_name)
    print('    ch2thf:', pool[CH2THF], ' nadpp:', pool[NADPP], ' chpthf:', pool[CHPTHF], ' nadph:', pool[NADPH])

    A = pool[CH2THF]
    B = pool[NADPP]
    
    #Constants for enzyme kinetics
    Vmax_2a = 1.0
    Km_ch2thf = 0.05
    Km_nadpp = 0.05

    # Calculate rate using two substrate MM equation
    rate_2a = two_substrate_mm(A, B, Vmax_2a, Km_ch2thf, Km_nadpp)
    delta_2a = rate_2a * step_timespan

    # New concentrations as reaction proceeds
    pool[CH2THF] -= delta_2a
    pool[NADPP] -= delta_2a
    pool[CHPTHF] += delta_2a
    pool[NADPH] += delta_2a
    
    print('    ch2thf:', pool[CH2THF], ' nadpp:', pool[NADPP], ' chpthf:', pool[CHPTHF], ' nadph:', pool[NADPH])
    return


def mito_chpthf_to_chothf(pool, step, step_timespan):
    """
    This function represents the second part of the second step in a biochemical pathway.
    It converts CH+THF into CHOTHF.

    Step 2b: [CH+THF] + [H2O] <--> [CHOTHF] + [H+]
    """
    print(' ', sys._getframe().f_code.co_name)
    print('    chpthf:', pool[CHPTHF], ' h2o:', pool[H2O], ' chothf:', pool[CHOTHF], ' h_plus:', pool[H_PLUS]) 
    
    Vmax_2b = 1.0
    Km_chpthf = 0.05
    Km_h2o = 0.05
    
    A = pool[CHPTHF]
    B = pool[H2O]
    
    # Calculate rate using two substrate MM equation
    rate_2b = two_substrate_mm(A, B, Vmax_2b, Km_chpthf, Km_h2o)
    delta_2b = rate_2b * step_timespan

    # New concentrations as reaction proceeds
    pool[CHPTHF] -= delta_2b
    pool[H2O] -= delta_2b
    pool[CHOTHF] += delta_2b
    pool[H_PLUS] += delta_2b

    print('    chpthf:', pool[CHPTHF], ' h2o:', pool[H2O], ' chothf:', pool[CHOTHF], ' h_plus:', pool[H_PLUS]) 
    return


def mito_chothf_atp_to_formate_adp(pool, step, step_timespan):
    """
    This function represents the third step in a biochemical pathway.
    It converts CHOTHF and ATP into Formate and ADP.

    Step 3: [CHOTHF] + [ATP] <--> [Formate] + [ADP]
    """
    print(' ', sys._getframe().f_code.co_name)
    print('    chothf:', pool[CHOTHF], ' atp:', pool[ATP], ' formate:', pool[FORMATE], ' adp:', pool[ADP])

    Vmax_3 = 1.0
    Km_chothf = 0.05
    Km_atp = 0.05
    
    A = pool[CHOTHF]
    B = pool[ATP]

    # Calculate rate using two substrate MM equation
    rate_3 = two_substrate_mm(A, B, Vmax_3, Km_chothf, Km_atp)
    delta_3 = rate_3 * step_timespan

    # New concentrations as reaction proceeds
    pool[CHOTHF] -= delta_3
    pool[ATP] -= delta_3
    pool[FORMATE] += delta_3
    pool[ADP] += delta_3
    
    print('    chothf:', pool[CHOTHF], ' atp:', pool[ATP], ' formate:', pool[FORMATE], ' adp:', pool[ADP])
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
        H_PLUS : 0.0,
        H2O : 0.1,
        CH2THF : 0.0,
        GLY : 0.0,
        FORMATE : 0.0,
        ADP : 0.0,
        CHOTHF : 0.0,
        CHPTHF : 0.0,
        NADPH : 0.0
    }

    print('Initial concentrations:')
    pprint(mito_pool)
    print()

    for step in range(0, experiment_timesteps):
        print('Step:', step, ' Timespan:', step_timespan)
        # Step 1: [SER] + [THF] <--> [CH2THF] + [GLY]
        mito_ser_thf_to_ch2thf_gly(mito_pool, step, step_timespan)

        # Step 2a: [CH2THF] + [NADP+] <--> [CH+THF] + [NADPH]
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
