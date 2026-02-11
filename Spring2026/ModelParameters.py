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

EXPERIMENT_RUNTIME_IN_SEC = 1  # Experiment time in seconds
EXPERIMENT_TIMESTEPS = 20  # Resolution for the simulation

INITIAL_SER_CONCENTRATION = 1963  # Concentration of SER in uM

# Initialize a pool of compounds with their concentrations
MITO_POOL = {
    SER: INITIAL_SER_CONCENTRATION,
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

