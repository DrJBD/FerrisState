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
NH3 = 'NH3'
CO2 = 'CO2'
H_OX = 'H_OX'
H_AM = 'H_AM'
H_RED = 'H_RED' # reduced H-protein: H-lip-(SH)2
NAD_PLUS = 'NAD+'
NADH = 'NADH'

STATE_COMPOUNDS = (
    SER,
    THF,
    NADPP,
    ATP,
    H_PLUS,
    H2O,
    CH2THF,
    GLY,
    FORMATE,
    ADP,
    CHOTHF,
    CHPTHF,
    NADPH,
    CO2,
    NH3,
    H_OX,
    H_AM,
    H_RED,
    NAD_PLUS,
    NADH
)

#----------Step 1--------------------------------------------------------------------------------------------------------------------------------------------------
# SHMT2 Reaction Parameters:
Step1_Vmax_f = 9756  # µM/h (Calculated. See word document)
Step1_Vmax_r = 3608  # µM/h (Calculated. See word document)
Step1_aKm_thf = 20  # µM (+/- 8)
Step1_Ki_thf = 99  # µM (+/- 17)
Step1_Ki_ser = 602  # µM (+/- 370)
Step1_aKm_gly = 469.8  # µM (Calculated. See word document)
Step1_Km_ch2thf = 980  # µM (Calculated. See word document)
Step1_KM_SER_MAX = 430.0      # µM
Step1_ALPHA_KM_SER = 158.0   # µM
Step1_K05_THF = 250.0        # µM

#---------GCS Sequential Reaction Parameters------------------------------------------------------------------------------------------------------------------------------------------------------
# Glycine Cleavage System:
# Overall: GLY + THF + NAD+ --> CH2THF + CO2 + NH3 + NADH + H+

# P-protein:
# GLY + H_OX --> H_AM + CO2
GCS_P_Vmax = 13200.0      # µM/h,
GCS_P_KM_GLY = 6600.0     # µM
GCS_P_KM_HOX = 7.4        # µM

# T-protein:
# H_AM + THF --> H_RED + CH2THF + NH3
GCS_T_Vmax = 44600.0      # µM/h,
GCS_T_KM_HAM = 584        # µM
GCS_T_KM_THF = 182100     # µM

# L-protein:
# H_RED + NAD+ --> H_OX + NADH + H+
GCS_L_Vmax = 13041     # µM/h
GCS_L_KM_HRED = 490    # µM
GCS_L_KM_NAD = 490     # µM


#------------Step 2a------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# MTHFD2 and MTHFD2L are treated as independent enzyme populations.
#
# Within each enzyme population, NAD+ and NADP+ are treated as
# alternative competing cofactors.
#
# Experimental mixed-cofactor inhibition constants (Ki) are not
# currently available.
#
# Therefore, as a baseline modeling assumption:
#
# Ki_competing_cofactor ≈ Km_competing_cofactor
#
# IMPORTANT:
# Km and Ki are NOT assumed to be inherently equivalent.
# Separate Ki variables are defined so experimentally measured Ki values can replace these assumed values
# later without changing the kinetic equations.
# ------------------------------------------------------------------
MTHFD2_NAD_VMAX = 477200.0           # µM/h
MTHFD2_NAD_KM_CH2THF = 354.0         # µM
MTHFD2_NAD_KM_NAD = 202.0            # µM
MTHFD2_NAD_KI_NADP = 352.0           # µM, assumed

MTHFD2_NADP_VMAX = 1152000.0         # µM/h
MTHFD2_NADP_KM_CH2THF = 302.0        # µM
MTHFD2_NADP_KM_NADP = 352.0          # µM
MTHFD2_NADP_KI_NAD = 202.0           # µM, assumed

# ==================================================================
MTHFD2L_NAD_VMAX = 1584000.0         # µM/h
MTHFD2L_NAD_KM_CH2THF = 130.0        # µM
MTHFD2L_NAD_KM_NAD = 147.0           # µM
MTHFD2L_NAD_KI_NADP = 537.0          # µM, assumed

MTHFD2L_NADP_VMAX = 1296000.0        # µM/h
MTHFD2L_NADP_KM_CH2THF = 153.0       # µM
MTHFD2L_NADP_KM_NADP = 537.0         # µM
MTHFD2L_NADP_KI_NAD = 147.0          # µM, assumed

#-------------Step 2b------------------------------------------------------------------------------------------------------------------------------------------------------------
#MTHFD2 Reaction Parameters
Step2b_Vmax_f = 8000  # µM/h
Step2b_Km_chpthf = 60  # µM
Step2b_Km_h2o = 40000  # µM (approx. cellular water conc.)

#-------------Step 3------------------------------------------------------------------------------------------------------------------------------------------------------------
#MTHFD1L Reaction Parameters
# Constants reverse
Step3_Vmax_f = 12000  # µM/h
Step3_Km_chothf = 100  # µM
Step3_Km_atp = 300  # µM
Step3_Km_adp = 300  #Temporary ?--------------------

#-------------Step 4------------------------------------------------------------------------------------------------------------------------------------------------------------
#ALDH1L2 Reaction Parameters
# Constants reverse
Step4_Vmax = 8640.0  # µM / hour---------------------To be changed (29.5 nmol/min/mg)
Step4_Km_chothf = 11.4  # µM


EXPERIMENT_RUNTIME_IN_HOURS = 1  # Experiment time will be divided by the number of timesteps
EXPERIMENT_TIMESTEPS = 1800   # Number of samples in the experiment (e.g., 1800 for 1 sample every 2 seconds over 1 hour)

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
    CHOTHF: 50, #doublecheck literature
    CHPTHF: 1.407,
    NADPH: 200,
    CO2: 0.0, #Need a literature value
    NH3: 0.0, #Need a literature value
    H_OX: 10.0, #Need a literature value
    H_AM: 0.0,
    H_RED: 0.0,
    NAD_PLUS: 100.0, #find ratio in literature
    NADH: 0.0 #Find ration in literature
}

