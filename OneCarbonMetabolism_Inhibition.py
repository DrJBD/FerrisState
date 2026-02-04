"""This script serves as a template for a Python program that includes a main function."""
import sys
from pprint import pprint
import matplotlib.pyplot as plt

# Reports
net_rate1_ch2thf = []
net_rate2_ch2thf = []
net_rate_gly = []
net_rate_report = []
net_rate_formate = []
net_rate_chothf = []
net_rate_chpthf = []


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
    delta = net_rate * step_timespan

    # Update pool
    pool[SER] -= delta
    pool[THF] -= delta
    pool[CH2THF] += delta
    pool[GLY] += delta

    # Clamp SER and THF after update
    pool[SER] = 100
    pool[THF] = 45
    pool[CH2THF] += 25  # µM per timestep
    pool[GLY] -= 40     # µM per timestep
    pool[NADPP] = 10000

    net_rate1_ch2thf.append(pool[CH2THF])
    net_rate_gly.append(pool[GLY])
    net_rate_report.append(net_rate)

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

    net_rate2_ch2thf.append(pool[CH2THF])

    #Clamping the metabolites in excess.

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

    net_rate_chpthf.append(pool[CHPTHF])
    net_rate_chothf.append(pool[CHOTHF])

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

    net_rate_formate.append(pool[FORMATE])

    print('    chothf:', pool[CHOTHF], ' atp:', pool[ATP], ' formate:', pool[FORMATE], ' adp:', pool[ADP])


    return

def ser_to_formate_mm_equation(ser_conc, thf_conc):
    """
            This function calculates the concentration of Formate produced from a given concentration of SER

            ser_conc: Concentration of SER as input
            Returns the concentration of Formate produced.
            """
    experiment_runtime = 900 # Experiment time in seconds
    experiment_timesteps = 450  # Resolution for the simulation. Each timestep lasts 2 sec.
    step_timespan = experiment_runtime / experiment_timesteps

    # Initialize a pool of compounds with their concentrations
    mito_pool = {
        SER: ser_conc,
        THF: thf_conc,     #µM [5]
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

def generate_plots(x_items, y_items, plot_label, y_axis_label, title):
    """
    This function generates plots for the simulation results.
    """

    plt.figure()
    #for clamp, style in thf_conditions:
    #    #@ t_s, ch2 = simulate_ch2thf_with_thf_clamp_zigzag(clamp_thf=clamp)
    #    t_s =
    #    # keep all curves black but distinguish by pattern/markers
    #    kwargs = dict(color="black", linewidth=1.2, zorder=2)
    #    if style["marker"]:
    #        kwargs.update(marker=style["marker"], markevery=60, markersize=3)
    #    plt.plot(t_s, ch2, linestyle=style["linestyle"], label=style["label"], **kwargs)

    kwargs = dict(color="#CF9DCD", linewidth=1.2, zorder=2)
    kwargs.update(marker="o", markevery=60, markersize=3)
    plt.plot(x_items, y_items, linestyle="-", label=plot_label, **kwargs,)

    # reference line at zero to highlight negative CH2THF if/when it occurs
    plt.axhline(0, linewidth=0.8, color="black", zorder=1)

    plt.xlabel("Time (s)")
    plt.ylabel(plot_label)
    plt.title(title)
    plt.legend()
    plt.tight_layout()
    plt.show()


def main():
    """
    This is the main function that will be executed when the script is run.
    """

    thf_concentrations = [20.8, 45, 100, 200]
    ser_concentration = 1963 #µM

    for thf_concentration in thf_concentrations:
        f = ser_to_formate_mm_equation(ser_concentration, thf_concentration)
        print('Formate produced from', ser_concentration, 'M SER:', f, 'M')

        print('net_rate_ch2thf 1:', net_rate1_ch2thf)
        print('net_rate_ch2thf 2:', net_rate2_ch2thf)
        print('net_rate_gly:', net_rate_gly)
        print('net_rate:', net_rate_report)
        print('net_rate_formate:', net_rate_formate)
        print('net_rate_chothf:', net_rate_chothf)
        print('net_rate_chpthf:', net_rate_chpthf)

        # Plotting the results
        points = 900
        x_steps = list(range(points))
        generate_plots(x_steps,
                       net_rate1_ch2thf[:points],
                       'net rate1 ch2thf ' + str(thf_concentration),
                       'net_rate1_ch2thf',
                       'Effect of [SER]= 100 µM on CH2THF production')
        generate_plots(x_steps, net_rate2_ch2thf[:points],
                       'net rate2 ch2thf ' + str(thf_concentration),
                       'net_rate2_ch2thf',
                       'Effect of [THF] TWO on CH2THF')
        generate_plots(x_steps, net_rate_gly[:points],
                       'net_rate_gly ' + str(thf_concentration),
                       'net_rate_gly',
                       'Effect of [[SER]= 2500 µM on Glycine')
        generate_plots(x_steps, net_rate_report[:points],
                       'net_rate_report ' + str(thf_concentration),
                       'net_rate_report',
                       'Effect of [SER]= 2500 µM on Net Reaction Rate')
        generate_plots(x_steps, net_rate_formate[:points],
                       'net_rate_formate ' + str(thf_concentration),
                          'net_rate_formate',
                       'Formate Production Over Time')
        generate_plots(x_steps, net_rate_chothf[:points],
                       'net_rate_chothf ' + str(thf_concentration),
                       'net_rate_chothf',
                       'Effect of [SER]= 2500 µM on CHOTHF production')
        generate_plots(x_steps, net_rate_chpthf[:points],
                       'net_rate_chpthf ' + str(thf_concentration),
                          'net_rate_chpthf',
                       'Effect of [SER]= 2500 µM on CH+THF production')

# Execute the main function
main()
