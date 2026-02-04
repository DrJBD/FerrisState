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
CO2= 'CO2'

#Constants
KM_SER_MAX = 430.0      # µM
ALPHA_KM_SER = 158.0   # µM
K05_THF = 250.0        # µM

def km_ser_app_from_thf(thf):
    """
    Apparent Km for serine as a function of THF concentration.
    Based on Fig. 2F.

    thf: [THF] in µM
    returns: Km_ser_app in µM
    """
    return (
        ALPHA_KM_SER
        + (KM_SER_MAX - ALPHA_KM_SER) * K05_THF / (K05_THF + thf)
    )

def shmt2_forward_inhibited(thf, ser, Vmaxf, aKm_thf, Ki_thf, Km_ser_app, Ki_ser):
    """
    Calculates the rate of the SHMT2 forward reaction with substrate inhibition.
    Ser + THF --> CH2THF + Gly
    """
    thf_term = thf / (aKm_thf + thf * (1 + thf / Ki_thf))
    ser_term = ser / (Km_ser_app + ser * (1 + ser / Ki_ser))
    return Vmaxf * thf_term * ser_term


def shmt2_reverse_rate(gly, ch2thf, Vmaxr, aKm_gly, Km_ch2thf):
    """
    Calculates the rate of the SHMT2 reverse reaction with Michaelis-Menten kinetics.
    Gly + CH2THF --> Ser + THF
    """
    gly_term = gly / (aKm_gly + gly)
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
    Vmax_f = 9756     #µM/h (Calculated. See word document)
    Vmax_r = 3608     #µM/h (Calculated. See word document)
    aKm_thf = 20      #µM (+/- 8)
    Ki_thf = 99       #µM (+/- 17)
    Ki_ser = 602      #µM (+/- 370)
    aKm_gly = 469.8   #µM (Calculated. See word document)
    Km_ch2thf = 980   #µM (Calculated. See word document)

    # Concentrations
    thf = pool[THF]
    Km_ser_app = km_ser_app_from_thf(thf)
    ser = pool[SER]
    gly = pool[GLY]
    ch2thf = pool[CH2THF]

    # Rates
    rate_forward = shmt2_forward_inhibited(thf, ser, Vmax_f, aKm_thf, Ki_thf, Km_ser_app, Ki_ser)
    rate_reverse = shmt2_reverse_rate(gly, ch2thf, Vmax_r, aKm_gly, Km_ch2thf)
    net_rate = rate_forward - rate_reverse
    dt_hr = step_timespan / 3600.0
    delta = net_rate * dt_hr

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
    Calculates the net velocity of a reversible (F&R) two-substrate reaction (S&F).
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

    # Constants forward
    Vmax_f = 2106 #µM/h [5]
    Km_ch2thf = 153  #µM [6] Pentaglutamated
    Km_nadpp = 537  #µM [6]

    #Constants reverse
    Vmax_r = 8000
    Km_chpthf = 80
    Km_nadph = 537 #µM

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

    dt_hr = step_timespan / 3600.0
    delta = net_rate * dt_hr

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
    Vmax_f = 8000  # µM/h
    Km_chpthf = 60  # µM
    Km_h2o = 40000  # µM (approx. cellular water conc.)

    # Constants reverse
    Vmax_r = 5000  # µM/h
    Km_chothf = 120  # µM
    Km_hplus = 10000  # µM (estimate)

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

    dt_hr = step_timespan / 3600.0
    delta = net_rate * dt_hr

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
        Reaction 3: CHOTHF + ATP <-->  Formate + THF + ADP
        Modeled using reversible two-substrate Michaelis-Menten kinetics.
        """

    print(' ', sys._getframe().f_code.co_name)
    print('    chothf:', pool[CHOTHF], ' atp:', pool[ATP], ' formate:', pool[FORMATE],' adp:', pool[ADP])

    # Constants (approximate; adjust based on literature)
    Vmax_f = 12000  # µM/h
    Km_chothf = 100  # µM
    Km_atp = 300  # µM

    # Constants reverse
    Vmax_r = 7000  # µM/h
    Km_formate = 50  # µM
    Km_adp = 100  # µM

    # Concentrations
    chothf = pool[CHOTHF]
    atp = pool[ATP]
    formate = pool[FORMATE]
    thf = pool[THF]
    adp = pool[ADP]

    # Net rate
    net_rate = reversible_two_substrate_mm(
        S_f=chothf, F_f=atp, Vmax_f=Vmax_f, Km_Sf=Km_chothf, Km_Ff=Km_atp,
        S_r=formate, F_r=adp, Vmax_r=Vmax_r, Km_Sr=Km_formate, Km_Fr=Km_adp
    )

    dt_hr = step_timespan / 3600.0
    delta = net_rate * dt_hr

    # Update pool
    pool[CHOTHF] -= delta
    pool[ATP] -= delta
    pool[FORMATE] += delta
    pool[THF] += delta
    pool[ADP] += delta

    print('    rate_net:', net_rate)
    print('    chothf:', pool[CHOTHF], ' atp:', pool[ATP], ' formate:', pool[FORMATE], ' adp:', pool[ADP])
    return

def mito_fdh_dehydrogenase(pool, step, step_timespan):
    """
    FDH dehydrogenase activity:
    Formate + NADP+ -> CO2 + NADPH
    Irreversible two-substrate Michaelis-Menten
    """

    # Parameters
    Vmax = 8640.0       # µM / hour
    Km_formate = 3.2    # µM
    Km_NADP = 1.1       # µM

    # Substrate concentrations
    formate = pool[FORMATE]
    nadp = pool[NADPP]

    # Rate law
    rate = (
        Vmax
        * (formate / (Km_formate + formate))
        * (nadp / (Km_NADP + nadp))
    )

    dt_hr = step_timespan / 3600.0
    delta = rate * dt_hr

    # Update pools
    pool[FORMATE] -= delta
    pool[NADPP]    -= delta
    pool[CO2]     += delta
    pool[NADPH]   += delta

    return rate


def ser_to_formate_mm_equation(ser_conc):
    """
            This function calculates the concentration of Formate produced from a given concentration of SER

            ser_conc: Concentration of SER as input
            Returns the concentration of Formate produced.
            """
    experiment_runtime = 1800  # Experiment time in seconds
    experiment_timesteps = 900  # Resolution for the simulation
    step_timespan = experiment_runtime / experiment_timesteps

    # Initialize a pool of compounds with their concentrations
    mito_pool = {
        SER: ser_conc,
        THF: 20.8,     #µM [5]
        NADPP: 1000,   #NADP+:NADPH 1:10 ([6] says 1.0 mM NAD+ was used).
        ATP: 300,      #ATP:ADP 3:1
        H_PLUS: 10000, #needs to be in excess
        H2O: 40000,    #needs to be in excess
        CH2THF: 1.7,   #[5] other sources give a range from 2.5-25 µM
        GLY: 1858,     #µM [5]
        FORMATE: 0.0,  #Output
        ADP: 100,      #ATP:ADP 3:1
        CHOTHF: 16,    #µM [5]
        CHPTHF: 1.55,  #µM [5]
        NADPH: 10000,   #NADP+:NADPH 1:10
        CO2:0.0,

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

        # Step 3: [CHOTHF] + [ATP] <--> [Formate] + [ADP] + [THF]
        mito_chothf_atp_to_formate_adp(mito_pool, step, step_timespan)

        #Step 4: [CHOTHF] + [NADP+] <--> [CO2] + [NADPH] + [THF]
        mito_fdh_dehydrogenase(mito_pool, step, step_timespan)

        print()

    print('Final concentrations after the experiment:')
    pprint(mito_pool)
    return mito_pool[FORMATE]  # Return the concentration of Formate produced


def main():
    """
    This is the main function that will be executed when the script is run.
    """
    ser_concentration = 1963 #µM
    f = ser_to_formate_mm_equation(ser_concentration)
    print('Formate produced from', ser_concentration, 'M SER:', f, 'M')

    # Execute the main function

main()
