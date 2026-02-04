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
    Vmax_f = 9756  # µM/h
    Vmax_r = 3608  # µM/h
    aKm_thf = 20  # µM
    Ki_thf = 99  # µM
    Km_ser_app = 158  # µM (This is only true about high concentrations of the inhibitory molecule (THF). In low conc. the Km becomes  variable of the [THF] function
    Ki_ser = 602  # µM
    aKm_gly = 469.8  # µM
    Km_ch2thf = 980  # µM

    # Concentrations
    thf = pool[THF]
    ser = pool[SER]
    gly = pool[GLY]
    ch2thf = pool[CH2THF]

    # Rates
    rate_forward = shmt2_forward_inhibited(thf, ser, Vmax_f, aKm_thf, Ki_thf, Km_ser_app, Ki_ser)
    rate_reverse = shmt2_reverse_rate(gly, ch2thf, Vmax_r, aKm_gly, Km_ch2thf)
    net_rate = rate_forward - rate_reverse
    delta = net_rate * step_timespan

    # Update pool
    pool[SER] -= delta
    pool[THF] -= delta
    pool[CH2THF] += delta
    pool[GLY] += delta

    # Clamp SER and THF after update
    pool[SER] = 1963
    pool[THF] = 45
    pool[CH2THF] += 25  # µM per timestep
    pool[GLY] -= 40  # µM per timestep
    pool[NADPP] = 10000
    print('    rate_fwd:', rate_forward, ' rate_rev:', rate_reverse, ' net_rate:', net_rate)
    print('    ser:', pool[SER], ' thf:', pool[THF], ' ch2thf:', pool[CH2THF], ' gly:', pool[GLY])
    return


def two_substrate_mm(A, B, Km_A, Km_B, Vmax):
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

    A = pool[CH2THF]  # uM (This pool should resemble the updated concentrations (pool[CH2THF] += delta_1)
    B = pool[NADPP]  # uM Approximation
    # pool[CHPTHF]    #Initial intermediate concentration (uM)
    # pool[NADPH]      #uM Approximation

    # Constants for enzyme kinetics
    Vmax_2a = 15000  # unknown for the moment. Catalytic efficiency "kcat/Km" is 0.067 (1/sM). The given value is an approximation
    Km_ch2thf = 50  # uM
    Km_nadpp = 537  # uM

    # Calculate rate using two substrate MM equation
    rate_2a = two_substrate_mm(A, B, Vmax_2a, Km_ch2thf, Km_nadpp)
    delta_2a = rate_2a * step_timespan

    # New concentrations as reaction proceeds
    pool[CH2THF] -= delta_2a
    pool[NADPP] -= delta_2a
    pool[CHPTHF] += delta_2a
    pool[NADPH] += delta_2a

    # Clamping the metabolites in excess.

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

    Vmax_2b = 20000  # uM/h
    Km_chpthf = 50  # uM Approximation
    Km_h2o = 0.001  # uM Approximation

    A = pool[CHPTHF]  # updated values from "pool[CHPTHF] += delta_2a"
    B = pool[H2O]  # Excess (uM). Approximation
    # pool[CHOTHF]        #Initial concentration
    # pool[H_PLUS]        #Excess (uM). Approximation.

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

    Vmax_3 = 11440  # uM/h Approximation.
    Km_chothf = 50  # uM Approximation
    Km_atp = 40  # uM

    A = pool[CHOTHF]  # Updated (uM)
    B = pool[ATP]  # Excess (uM) Approximation
    # pool[FORMATE]          #(uM) If we're trying to isolate the net production of formate
    # pool[ADP]              #Excess (uM) Approximation

    # Calculate rate using two substrate MM equation
    rate_3 = two_substrate_mm(A, B, Vmax_3, Km_chothf, Km_atp)
    delta_3 = rate_3 * step_timespan

    # New concentrations as reaction proceeds
    pool[CHOTHF] -= delta_3
    pool[ATP] -= delta_3
    pool[FORMATE] += delta_3  # Outcome
    pool[ADP] += delta_3

    print('    chothf:', pool[CHOTHF], ' atp:', pool[ATP], ' formate:', pool[FORMATE], ' adp:', pool[ADP])
    return


def ser_to_formate_mm_equation(ser_conc):
    """
            This function calculates the concentration of Formate produced from a given concentration of SER

            ser_conc: Concentration of SER as input
            Returns the concentration of Formate produced.
            """
    experiment_runtime = 1800  # Experiment time in seconds
    experiment_timesteps = 900  # Resolution for the simulation. Each timestep lasts 2 sec.
    step_timespan = experiment_runtime / experiment_timesteps

    # Initialize a pool of compounds with their concentrations
    mito_pool = {
        SER: ser_conc,
        THF: 45,  # µM [5]
        NADPP: 10000,  # NADP+:NADPH 10:1
        ATP: 300,  # ATP:ADP 3:1
        H_PLUS: 10000,  # needs to be in excess
        H2O: 40000,  # needs to be in excess
        CH2THF: 25,  # [5] other sources give a range from 2.5-25 µM
        GLY: 1.857,  # µM [5]
        FORMATE: 0.0,  # Output
        ADP: 100,  # ATP:ADP 3:1
        CHOTHF: 16,  # µM [5]
        CHPTHF: 1.55,  # µM [5]
        NADPH: 9500,  # NADP+:NADPH 1:10
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
    ser_concentration = 1963  # µM
    f = ser_to_formate_mm_equation(ser_concentration)
    print('Formate produced from', ser_concentration, 'M SER:', f, 'M')

    # Execute the main function


main()


# ===================== APPEND BELOW YOUR EXISTING main() CALL =====================


def simulate_ch2thf_with_thf_clamp_zigzag(clamp_thf, ser_conc=1963,
                                          experiment_runtime=1800, experiment_timesteps=900):
    """
    Run one simulation while forcing THF to `clamp_thf` at every step.
    We record CH2THF twice per step to capture the zig-zag:
      - after Step 1 (mid-step time)
      - after Step 3 (end-of-step time)
    Returns (times_seconds, ch2thf_series_uM).
    """
    step_timespan = experiment_runtime / experiment_timesteps  # seconds

    # Initial pool (same as your function, but THF starts at the clamp)
    pool = {
        SER: ser_conc,
        THF: float(clamp_thf),
        NADPP: 10000,
        ATP: 300,
        H_PLUS: 10000,
        H2O: 40000,
        CH2THF: 25,
        GLY: 1.857,
        FORMATE: 0.0,
        ADP: 100,
        CHOTHF: 16,
        CHPTHF: 1.55,
        NADPH: 9500,
    }

    times = [0.0]
    ch2_series = [pool[CH2THF]]
    t = 0.0

    for step in range(experiment_timesteps):
        # --- Step 1 ---
        mito_ser_thf_to_ch2thf_gly(pool, step, step_timespan)

        # Immediately re-clamp THF for this experiment
        pool[THF] = float(clamp_thf)

        # Record mid-step (after Step 1) to capture the upward spike
        times.append(t + 0.5 * step_timespan)
        ch2_series.append(pool[CH2THF])

        # --- Steps 2a, 2b, 3 ---
        mito_ch2thf_nadpp_to_chpthf_nadph(pool, step, step_timespan)
        mito_chpthf_to_chothf(pool, step, step_timespan)
        mito_chothf_atp_to_formate_adp(pool, step, step_timespan)

        # End of this timestep
        t += step_timespan
        times.append(t)
        ch2_series.append(pool[CH2THF])

        # Ensure clamp is in place before the next loop iteration
        pool[THF] = float(clamp_thf)

    return times, ch2_series


# ---- Plot CH2THF vs time for four THF clamps (all black, different patterns) ----
thf_conditions = [
    (20.8, {"linestyle": "-", "marker": None, "label": "THF=20.8 µM"}),  # solid
    (45, {"linestyle": ":", "marker": None, "label": "THF=45 µM"}),  # dotted
    (100.0, {"linestyle": "-", "marker": "^", "label": "THF=100 µM"}),  # triangles
    (200.0, {"linestyle": "-", "marker": "o", "label": "THF=200 µM"}),  # circles
]

plt.figure()
for clamp, style in thf_conditions:
    t_s, ch2 = simulate_ch2thf_with_thf_clamp_zigzag(clamp_thf=clamp)
    # keep all curves black but distinguish by pattern/markers
    kwargs = dict(color="black", linewidth=1.2, zorder=2)
    if style["marker"]:
        kwargs.update(marker=style["marker"], markevery=60, markersize=3)
    plt.plot(t_s, ch2, linestyle=style["linestyle"], label=style["label"], **kwargs)

# reference line at zero to highlight negative CH2THF if/when it occurs
plt.axhline(0, linewidth=0.8, color="black", zorder=1)

plt.xlabel("Time (s)")
plt.ylabel("CH2THF (µM)")
plt.title("CH2THF vs Time for Different THF Clamps (THF forced each step)")
plt.legend()
plt.tight_layout()
plt.show()
