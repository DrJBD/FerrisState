def glycine_cleavage_system(pool, step_timespan):
    """
    Simulated GCS reaction (one-step approximation (maybe)):
    Glycine + THF + NAD+ → CH2THF + CO2 + NH3 + NADH + H+
    Using placeholder kinetics—adjust values when real data is available.
    """
    gly = pool[GLY]
    thf = pool[THF]
    nadp+ = pool[NADPP]

    # Placeholder kinetic constants. To be expanded
    Vmax = 400      # adjust later (units: µM per time unit)
    Km_gly = 500    # placeholder (µM)
    Km_thf = 600     # placeholder (µM)
    Km_nad = 700     # placeholder (µM)

    # Rate Law to be determined.
    #rate = Vmax * (gly / (Km_gly + gly)) * (thf / (Km_thf + thf)) * (nadp / (Km_nad + nadp))
    delta = rate * step_timespan

    # Update pools
    pool[GLY]   -= delta
    pool[THF]   -= delta
    pool[NADPP] -= delta
    pool[CH2THF] += delta
    pool[NADPH] += delta
    pool[H_PLUS] = pool.get(H_PLUS, 0) + delta
    pool['CO2']  = pool.get('CO2', 0) + delta
    pool['NH3']  = pool.get('NH3', 0) + delta

    return rate
