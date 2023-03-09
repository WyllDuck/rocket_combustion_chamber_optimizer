""" GLOBAL PARAMETERS """
# FIXED GLOBAL PARAMETERS
N_SECTIONS      = 300    # [º]       - number of cross-sections of the combustion chamber - simulation resolution

# COOLANT PARAMETERS
INLET_T         = 200   # [K]       - inlet temperature coolant
INLET_p         = 130   # [bar]     - inlet pressure coolant
MDOT_COOLANT    = 80    # [kg/s]    - mass flow rate of the coolant

# MATERIAL PARAMETERS
THERMAL_K       = 385   # [W/mK]    - thermal conductivity of the material
FRICTION        = 0.01  # [-]       - friction factor

# GEOMETRICAL PARAMETERS
N_CHANNELS      = 40    # [-]       - number of channels in the combustion chamber
INTER_CHANNEL_T = 0.005 # [m]       - thickness of the wall separating the channels
LENGHT_CC       = 1.0   # [m]       - length of the combustion chamber
HEIHT_CHANNEL   = 0.004 # [m]       - height of the cooling channel
DI              = 1.0   # [m]       - inner diameter of the combustion chamber
T               = 0.01  # [m]       - thickness of the wall separating the combustion chamber from the coolant circuit

# HOT GAS PARAMETERS
HOT_GAS_T       = 2000      # [K]       - temperature of the hot gases
HOT_GAS_CP      = 6950      # [J/kgK]   - specific heat of the hot gases
HOT_GAS_MU      = 1.17e-4   # [Pa.s]    - viscosity of the hot gases
HOT_GAS_K       = 1.49      # [W/mK]    - thermal conductivity of the hot gases
HOT_GAS_RHO     = 9.24      # [kg/m3]   - density of the hot gases
HOT_GAS_V       = 1237      # [m/s]     - velocity of the hot gases

THERMAL_PARA_HOT_GASES = [HOT_GAS_T, HOT_GAS_CP, HOT_GAS_MU, HOT_GAS_K, HOT_GAS_RHO, HOT_GAS_V]

# INITIAL ASSUMPTIONS
T_WI_INIT_ASSUMPTION = HOT_GAS_T - 1000 # [K] - initial temperature assumption for inner combustion chamber wall
T_WO_INIT_ASSUMPTION = HOT_GAS_T - 1200 # [K] - initial temperature assumption for outer combustion chamber wall

# COOLANT PROPERTIES
CONF_COOLANT_FILES = "conf/coolant"
