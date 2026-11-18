"""This script serves as a template for a Python program that includes a main function."""
import sys
from pprint import pprint
import matplotlib.pyplot as plt


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
    Vmax_f = 9756     #µM/h
    Vmax_r = 3608     #µM/h
    aKm_thf = 20      #µM
    Ki_thf = 99       #µM
    Km_ser_app = 158  #µM (This is only true about high concentrations of the inhibitory molecule (THF). In low conc. the Km becomes  variable of the [THF] function
    Ki_ser = 602      #µM
    aKm_gly = 469.8   #µM
    Km_ch2thf = 980   #µM

    # Concentrations
    thf = pool[THF]
    ser = pool[SER]
    gly = pool[GLY]
    ch2thf = pool[CH2THF]

    # Rates
    rate_forward = shmt2_forward_inhibited(thf, ser, Vmax_f, aKm_thf, Ki_thf, Km_ser_app, Ki_ser)
    rate_reverse = shmt2_reverse_rate(gly, ch2thf, Vmax_r, aKm_gly, Km_ch2thf)
    net_rate = rate_forward - rate_reverse
    delta = net_rate * step_timespan  # (keep your original scaling)

    # Update pool
    pool[SER] -= delta
    pool[THF] -= delta
    pool[CH2THF] += delta
    pool[GLY] += delta

    # Clamp SER and THF after update  (keep exactly as-is)
    pool[SER] = 100
    pool[THF] = 45
    pool[CH2THF] += 25  # µM per timestep
    pool[GLY] -= 40     # µM per timestep
    pool[NADPP] = 10000

    print('    rate_fwd:', rate_forward, ' rate_rev:', rate_reverse, ' net_rate:', net_rate)
    print('    ser:', pool[SER], ' thf:', pool[THF], ' ch2thf:', pool[CH2THF], ' gly:', pool[GLY])
    return


def two_substrate_mm(A, B, Vmax, Km_A, Km_B):
    """
    This function computes the rate for a two-substrate reaction
    A: Concentration of substrate A
    B: Concentration of substrate B
    Vmax: Maximum rate of reaction.
    Km_A, Km_B: Taken from literature. Constants. Dictate the concentration of substrates as half-Vmax is reached.
    """
    return (Vmax * A * B) / (Km_A * B + Km_B * A + A * B)


def mito_ch2thf_nadpp_to_chpthf_nadph(pool, step, step_timespan):
    """
    This function represents the second step in a biochemical pathway.
    It converts CH2THF and NADP+ into CH+THF and NADPH.

    Step 2a: [CH2THF] + [NADP+] <--> [CH+THF] + [NADPH]
    """

    print(' ', sys._getframe().f_code.co_name)
    print('    ch2thf:', pool[CH2THF], ' nadpp:', pool[NADPP], ' chpthf:', pool[CHPTHF], ' nadph:', pool[NADPH])

    A = pool[CH2THF]  # uM
    B = pool[NADPP]   # uM

    # Constants for enzyme kinetics
    Vmax_2a = 15000
    Km_ch2thf = 50
    Km_nadpp = 537

    # Keep your original call order
    rate_2a = two_substrate_mm(A, B, Vmax_2a, Km_ch2thf, Km_nadpp)
    delta_2a = rate_2a * step_timespan

    pool[CH2THF] -= delta_2a
    pool[NADPP]  -= delta_2a
    pool[CHPTHF] += delta_2a
    pool[NADPH]  += delta_2a

    print('    ch2thf:', pool[CH2THF], ' nadpp:', pool[NADPP], ' chpthf:', pool[CHPTHF], ' nadph:', pool[NADPH])
    return


def mito_chpthf_to_chothf(pool, step, step_timespan):
    """
    Step 2b: [CH+THF] + [H2O] <--> [CHOTHF] + [H+]
    """
    print(' ', sys._getframe().f_code.co_name)
    print('    chpthf:', pool[CHPTHF], ' h2o:', pool[H2O], ' chothf:', pool[CHOTHF], ' h_plus:', pool[H_PLUS])

    Vmax_2b = 20000
    Km_chpthf = 50
    Km_h2o = 0.001

    A = pool[CHPTHF]
    B = pool[H2O]

    # Keep your original call order
    rate_2b = two_substrate_mm(A, B, Vmax_2b, Km_chpthf, Km_h2o)
    delta_2b = rate_2b * step_timespan

    pool[CHPTHF] -= delta_2b
    pool[H2O]    -= delta_2b
    pool[CHOTHF] += delta_2b
    pool[H_PLUS] += delta_2b

    print('    chpthf:', pool[CHPTHF], ' h2o:', pool[H2O], ' chothf:', pool[CHOTHF], ' h_plus:', pool[H_PLUS])
    return


def mito_chothf_atp_to_formate_adp(pool, step, step_timespan):
    """
    Step 3: [CHOTHF] + [ATP] <--> [Formate] + [ADP]
    """
    print(' ', sys._getframe().f_code.co_name)
    print('    chothf:', pool[CHOTHF], ' atp:', pool[ATP], ' formate:', pool[FORMATE], ' adp:', pool[ADP])

    Vmax_3 = 11440
    Km_chothf = 50
    Km_atp = 40

    A = pool[CHOTHF]
    B = pool[ATP]

    # Keep your original call order
    rate_3 = two_substrate_mm(A, B, Vmax_3, Km_chothf, Km_atp)
    delta_3 = rate_3 * step_timespan

    pool[CHOTHF] -= delta_3
    pool[ATP]    -= delta_3
    pool[FORMATE] += delta_3
    pool[ADP]    += delta_3

    print('    chothf:', pool[CHOTHF], ' atp:', pool[ATP], ' formate:', pool[FORMATE], ' adp:', pool[ADP])
    return


# ---------- SHMT2 "snapshot" for plotting (does not mutate pool) ----------
def shmt2_rates_snapshot(pool):
    """Compute SHMT2 forward/reverse/net rates from the *current* pool without updating it."""
    Vmax_f = 9756
    Vmax_r = 3608
    aKm_thf = 20
    Ki_thf = 99
    Km_ser_app = 158
    Ki_ser = 602
    aKm_gly = 469.8
    Km_ch2thf = 980

    thf = pool[THF]
    ser = pool[SER]
    gly = pool[GLY]
    ch2thf = pool[CH2THF]

    rate_forward = shmt2_forward_inhibited(thf, ser, Vmax_f, aKm_thf, Ki_thf, Km_ser_app, Ki_ser)
    rate_reverse = shmt2_reverse_rate(gly, ch2thf, Vmax_r, aKm_gly, Km_ch2thf)
    net_rate = rate_forward - rate_reverse
    return rate_forward, rate_reverse, net_rate
# ------------------------------------------------------------------------


def ser_to_formate_mm_equation(ser_conc):
    """
            This function calculates the concentration of Formate produced from a given concentration of SER

            ser_conc: Concentration of SER as input
            Returns the concentration of Formate produced.
            """
    experiment_runtime = 1800 # Experiment time in seconds
    experiment_timesteps = 900  # Resolution for the simulation. Each timestep lasts 2 sec.
    step_timespan = experiment_runtime / experiment_timesteps

    # Initialize a pool of compounds with their concentrations
    mito_pool = {
        SER: ser_conc,
        THF: 45,     #µM [5]
        NADPP: 10000,   #NADP+:NADPH 10:1
        ATP: 300,      #ATP:ADP 3:1
        H_PLUS: 10000, #needs to be in excess
        H2O: 40000,    #needs to be in excess
        CH2THF: 25,    #[5] other sources give a range from 2.5-25 µM
        GLY: 1.857,    #µM [5]
        FORMATE: 0.0,  #Output
        ADP: 100,      #ATP:ADP 3:1
        CHOTHF: 16,    #µM [5]
        CHPTHF: 1.55,  #µM [5]
        NADPH: 9500,    #NADP+:NADPH 1:10
    }

    # Recorders for plotting (do not influence dynamics)
    times = [0.0]                         # seconds
    formate_series = [mito_pool[FORMATE]] # µM
    shmt2_fwd_rates = []                  # µM/h per timestep
    shmt2_rev_rates = []                  # µM/h per timestep
    shmt2_net_rates = []                  # µM/h per timestep

    print('Initial concentrations:')
    pprint(mito_pool)
    print()

    for step in range(0, experiment_timesteps):
        print('Step:', step, ' Timespan:', step_timespan)

        # Snapshot BEFORE mutating Step 1
        r_fwd, r_rev, r_net = shmt2_rates_snapshot(mito_pool)
        shmt2_fwd_rates.append(r_fwd)
        shmt2_rev_rates.append(r_rev)
        shmt2_net_rates.append(r_net)

        # Step 1
        mito_ser_thf_to_ch2thf_gly(mito_pool, step, step_timespan)
        # Step 2a
        mito_ch2thf_nadpp_to_chpthf_nadph(mito_pool, step, step_timespan)
        # Step 2b
        mito_chpthf_to_chothf(mito_pool, step, step_timespan)
        # Step 3
        mito_chothf_atp_to_formate_adp(mito_pool, step, step_timespan)
        print()

        # Record for the formate time-course
        current_time = (step + 1) * step_timespan  # seconds
        times.append(current_time)
        formate_series.append(mito_pool[FORMATE])

    print('Final concentrations after the experiment:')
    pprint(mito_pool)
    return (
        mito_pool[FORMATE],  # final formate
        times, formate_series,
        shmt2_fwd_rates, shmt2_rev_rates, shmt2_net_rates
    )


def main():
    """
    This is the main function that will be executed when the script is run.
    """
    ser_concentration =3000 #µM
    (final_formate,
     times, formate_series,
     shmt2_fwd_rates, shmt2_rev_rates, shmt2_net_rates) = ser_to_formate_mm_equation(ser_concentration)

    print('Formate produced from', ser_concentration, 'µM SER:', final_formate, 'µM')

    # ---- Plot 1: Formate over time ----
    plt.figure()
    plt.plot(times, formate_series, color="#CF9DCD")
    plt.xlabel('Time (s)')
    plt.ylabel('Formate (µM)')
    plt.title('Formate Production Over Time')
    plt.tight_layout()

    # ---- Plot 2: SHMT2 net-rate bars (first up to 70 steps) ----
    L = len(shmt2_net_rates)
    if L > 0:
        N = min(70, L)        # show up to ~70 timesteps
        x = list(range(1, N + 1))
        y = shmt2_net_rates[:N]
        width = 0.8

        plt.figure()
        # Black rectangles
        plt.bar(x, y, width=width, color='#CF9DCD')

        # Build right-edge connector lines, separately for positive and negative bars
        pos = [(xi, yi) for xi, yi in zip(x, y) if yi >= 0]
        neg = [(xi, yi) for xi, yi in zip(x, y) if yi < 0]

        if pos:
            pos.sort(key=lambda t: t[0])
            rx_pos = [xi + width/2 for xi, _ in pos]
            ry_pos = [yi for _, yi in pos]
            plt.plot(rx_pos, ry_pos, linewidth=1)  # line across tops of positive bars

        if neg:
            neg.sort(key=lambda t: t[0])
            rx_neg = [xi + width/2 for xi, _ in neg]
            ry_neg = [yi for _, yi in neg]
            plt.plot(rx_neg, ry_neg, linewidth=1)  # line across tops of negative (upside-down) bars

        plt.axhline(0, linewidth=1)
        plt.xlabel('Timestep')
        plt.ylabel('SHMT2 Net Rate (µM/h)')
        plt.title('SHMT2 Net Reaction Rate per Timestep when [SER]= 100')
        plt.xticks(x)
        plt.tight_layout()

        # ---- Plot 3: SHMT2 forward (+) vs reverse (−) (unchanged: first 25 steps) ----
        N_fr = min(25, L)
        x_fr = list(range(1, N_fr + 1))
        y_fwd = shmt2_fwd_rates[:N_fr]
        y_rev = [-v for v in shmt2_rev_rates[:N_fr]]  # show reverse as negative
        plt.figure()
        plt.bar(x_fr, y_fwd, label='Forward')
        plt.bar(x_fr, y_rev, label='Reverse')
        plt.axhline(0, linewidth=1)
        plt.xlabel('Timestep')
        plt.ylabel('Rate (µM/h)')
        plt.title('SHMT2 Forward (+) vs Reverse (−) Rates (First 25)')
        plt.xticks(x_fr)
        plt.legend()
        plt.tight_layout()
    else:
        print("No SHMT2 rates recorded — skipping bar charts.")

    plt.show()

# Execute the main function
main()
