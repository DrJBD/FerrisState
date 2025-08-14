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


def shmt2_forward_inhibited(thf, ser, Vmaxf, Km_thf, Ki_thf, Km_ser_app, Ki_ser, alpha=1):
    """
    Calculates the rate of the SHMT2 forward reaction with substrate inhibition.
    Ser + THF --> CH2THF + Gly
    """
    thf_term = thf / (alpha * Km_thf + thf * (1 + thf / Ki_thf))
    ser_term = ser / (Km_ser_app + ser * (1 + ser / Ki_ser))
    return Vmaxf * thf_term * ser_term

def shmt2_reverse_rate(gly, ch2thf, Vmaxr, Km_gly, Km_ch2thf, alpha=1):
    """
    Calculates the rate of the SHMT2 reverse reaction with Michaelis-Menten kinetics.
    Gly + CH2THF --> Ser + THF
    """
    gly_term = gly / (alpha * Km_gly + gly)
    ch2thf_term = ch2thf / (Km_ch2thf + ch2thf)
    return Vmaxr * gly_term * ch2thf_term

def mito_ser_thf_to_ch2thf_gly(pool, step, step_timespan):
        """
        This function models the reversible reaction:
        [SER] + [THF] --> [CH2THF] + [GLY]
        """

        print(' ', sys._getframe().f_code.co_name)
        print('    ser:', pool[SER], ' thf:', pool[THF], ' ch2thf:', pool[CH2THF], ' gly:', pool[GLY])

        # Constants
        Vmax_f = 11440  # µM/h
        Vmax_r = 8500  # µM/h (estimate)
        alpha = 1
        Km_thf = 50
        Ki_thf = 500
        Km_ser_app = 600
        Ki_ser = 2000
        Km_gly = 400
        Km_ch2thf = 100

        # Concentrations
        thf = pool[THF]
        ser = pool[SER]
        gly = pool[GLY]
        ch2thf = pool[CH2THF]

        # Rates
        rate_forward = shmt2_forward_inhibited(thf, ser, Vmax_f, Km_thf, Ki_thf, Km_ser_app, Ki_ser, alpha)
        rate_reverse = shmt2_reverse_rate(gly, ch2thf, Vmax_r, Km_gly, Km_ch2thf, alpha)
        net_rate = rate_forward - rate_reverse
        delta = net_rate * step_timespan

        # Update pool
        pool[SER] -= delta
        pool[THF] -= delta
        pool[CH2THF] += delta
        pool[GLY] += delta

        print('    rate_fwd:', rate_forward, ' rate_rev:', rate_reverse, ' net_rate:', net_rate)
        print('    ser:', pool[SER], ' thf:', pool[THF], ' ch2thf:', pool[CH2THF], ' gly:', pool[GLY])
        return

def reversible_two_substrate_mm(
    S_f, F_f, Vmax_f, Km_Sf, Km_Ff,
    S_r, F_r, Vmax_r, Km_Sr, Km_Fr
):

    """
    Calculates the net velocity of a reversible two-substrate reaction.
    This is the general form.

    Parameters:
    - S_f, F_f: concentrations of substrates in the forward reaction
    - Vmax_f, Km_Sf, Km_Ff: forward reaction constants
    - S_r, F_r: concentrations of substrates in the reverse reaction
    - Vmax_r, Km_Sr, Km_Fr: reverse reaction constants

    """
    forward = Vmax_f * (S_f / (Km_Sf + S_f)) * (F_f / (Km_Ff + F_f))
    reverse = Vmax_r * (S_r / (Km_Sr + S_r)) * (F_r / (Km_Fr + F_r))
    return forward - reverse

def mito_ch2thf_nadpp_to_chpthf_nadph(pool, step, step_timespan):
    """
    Reaction 2a: CH2THF + NADP⁺ <--> CH⁺THF + NADPH
    Modeled using reversible two-substrate Michaelis-Menten kinetics.
    """

    print(' ', sys._getframe().f_code.co_name)
    print('    ch2thf:', pool[CH2THF], ' nadpp:', pool[NADPP], ' chpthf:', pool[CHPTHF], ' nadph:', pool[NADPH])

    # Constants
    Vmax_f = 15000       # µM/h (estimate)
    Km_ch2thf = 50       # µM
    Km_nadpp = 537       # µM

    Vmax_r = 8000        # µM/h (example estimate)
    Km_chpthf = 80       # µM (estimate)
    Km_nadph = 200       # µM (estimate)

    # Concentrations
    ch2thf = pool[CH2THF]
    nadpp = pool[NADPP]
    chpthf = pool[CHPTHF]
    nadph = pool[NADPH]

    # Net reversible rate
    net_rate = reversible_two_substrate_mm(
        S_f=ch2thf, F_f=nadpp, Vmax_f=Vmax_f, Km_Sf=Km_ch2thf, Km_Ff=Km_nadpp,
        S_r=chpthf, F_r=nadph, Vmax_r=Vmax_r, Km_Sr=Km_chpthf, Km_Fr=Km_nadph
    )

    delta = net_rate * step_timespan

    # Update pool
    pool[CH2THF] -= delta
    pool[NADPP] -= delta
    pool[CHPTHF] += delta
    pool[NADPH] += delta

    print('    rate_net:', net_rate)
    print('    ch2thf:', pool[CH2THF], ' nadpp:', pool[NADPP], ' chpthf:', pool[CHPTHF], ' nadph:', pool[NADPH])
    return

def mito_chpthf_to_chothf(pool, step, step_timespan):
    """
    Reaction 2b: CH+THF + H2O <-->  CHOTHF + H+
    Modeled using reversible two-substrate Michaelis-Menten kinetics.
    """

    print(' ', sys._getframe().f_code.co_name)
    print('    chpthf:', pool[CHPTHF], ' h2o:', pool[H2O], ' chothf:', pool[CHOTHF], ' h+ :', pool[H_PLUS])

    # Constants (example estimates; refine with experimental values if available)
    Vmax_f = 8000         # µM/h
    Km_chpthf = 60        # µM
    Km_h2o = 40000        # µM (approx. cellular water conc.)

    Vmax_r = 5000         # µM/h
    Km_chothf = 120       # µM
    Km_hplus = 10000      # µM (estimate)

    # Concentrations
    chpthf = pool[CHPTHF]
    h2o = pool[H2O]
    chothf = pool[CHOTHF]
    hplus = pool[H_PLUS]

    # Net rate using reversible two-substrate MM
    net_rate = reversible_two_substrate_mm(
        S_f=chpthf, F_f=h2o, Vmax_f=Vmax_f, Km_Sf=Km_chpthf, Km_Ff=Km_h2o,
        S_r=chothf, F_r=hplus, Vmax_r=Vmax_r, Km_Sr=Km_chothf, Km_Fr=Km_hplus
    )

    delta = net_rate * step_timespan

    # Update the pool
    pool[CHPTHF] -= delta
    pool[H2O] -= delta
    pool[CHOTHF] += delta
    pool[H_PLUS] += delta

    print('    rate_net:', net_rate)
    print('    chpthf:', pool[CHPTHF], ' h2o:', pool[H2O], ' chothf:', pool[CHOTHF], ' h+ :', pool[H_PLUS])
    return

def mito_chothf_atp_to_formate_adp(pool, step, step_timespan):
        """
        Reaction 3: CHOTHF + ATP <-->  Formate + ADP
        Modeled using reversible two-substrate Michaelis-Menten kinetics.
        """

        print(' ', sys._getframe().f_code.co_name)
        print('    chothf:', pool[CHOTHF], ' atp:', pool[ATP], ' formate:', pool[FORMATE], ' adp:', pool[ADP])

        # Constants (approximate; adjust based on literature)
        Vmax_f = 12000  # µM/h (forward)
        Km_chothf = 100  # µM
        Km_atp = 300  # µM

        Vmax_r = 7000  # µM/h (reverse)
        Km_formate = 50  # µM
        Km_adp = 100  # µM

        # Concentrations
        chothf = pool[CHOTHF]
        atp = pool[ATP]
        formate = pool[FORMATE]
        adp = pool[ADP]

        # Net rate
        net_rate = reversible_two_substrate_mm(
            S_f=chothf, F_f=atp, Vmax_f=Vmax_f, Km_Sf=Km_chothf, Km_Ff=Km_atp,
            S_r=formate, F_r=adp, Vmax_r=Vmax_r, Km_Sr=Km_formate, Km_Fr=Km_adp
        )

        delta = net_rate * step_timespan

        # Update pool
        pool[CHOTHF] -= delta
        pool[ATP] -= delta
        pool[FORMATE] += delta
        pool[ADP] += delta

        print('    rate_net:', net_rate)
        print('    chothf:', pool[CHOTHF], ' atp:', pool[ATP], ' formate:', pool[FORMATE], ' adp:', pool[ADP])
        return

def ser_to_formate_mm_equation(ser_conc):
            """
            This function calculates the concentration of Formate produced from a given concentration of SER

            ser_conc: Concentration of SER as input
            Returns the concentration of Formate produced.
            """
            experiment_runtime = 1  # Experiment time in seconds
            experiment_timesteps = 20  # Resolution for the simulation
            step_timespan = experiment_runtime / experiment_timesteps

            # Initialize a pool of compounds with their concentrations
            mito_pool = {
                 SER: ser_conc,
                 THF: 20.8,
                 NADPP: 20,
                 ATP: 300,
                 H_PLUS: 10000,
                 H2O: 40000,
                 CH2THF: 50,
                 GLY: 1858,
                 FORMATE: 0.0,
                 ADP: 100,
                 CHOTHF: 1963,
                 CHPTHF: 1.407,
                 NADPH: 200
            }

            print('Initial concentrations:')
            pprint(mito_pool)
            print()

            for step in range(0, experiment_timesteps):
                print('Step:', step, ' Timespan:', step_timespan)
             # Step 1: [SER] + [THF] <--> [CH2THF] + [GLY]
                mito_ser_thf_to_ch2thf_gly(mito_pool, step, step_timespan)

             # Step 2a: [CH2THF] + [NADP+] <--> [CH+THF] + [NADPH]
                mito_ch2thf_nadpp_to_chpthf_nadph(mito_pool, step, step_timespan)

             # Step 2b: [CH+THF] <--> [CHOTHF]
                mito_chpthf_to_chothf(mito_pool, step, step_timespan)

             # Step 3: [CHOTHF] + [ATP] <--> [Formate] + [ADP]
                mito_chothf_atp_to_formate_adp(mito_pool, step, step_timespan)
                print()

            print('Final concentrations after the experiment:')
            pprint(mito_pool)
            return mito_pool[FORMATE]  # Return the concentration of Formate produced


def main():
    """
    This is the main function that will be executed when the script is run.
    """
    ser_concentration = 1963
    f = ser_to_formate_mm_equation(ser_concentration)
    print('Formate produced from', ser_concentration, 'M SER:', f, 'M')

    # Execute the main function


main()
