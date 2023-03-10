""" GLOBAL PARAMETERS """
N_SECTIONS      = 300       # [º]       - number of cross-sections of the combustion chamber - simulation resolution

# COOLANT PARAMETERS
INLET_T         = 200       # [K]       - inlet temperature coolant
INLET_p         = 125       # [bar]     - inlet pressure coolant
MDOT            = 324       # [kg/s]    - mass flow rate of hot gas
OF_RATIO        = 3.6       # [-]       - oxygen to fuel ratio

MDOT_COOLANT    = MDOT / (OF_RATIO + 1) # [kg/s] - mass flow rate of the coolant
"""
M_OXY / M_FUEL = OF_RATIO
MDOT = M_OXY + M_FUEL
MDOT = (OF_RATIO + 1) * M_FUEL
M_FUEL = MDOT_COOLANT / (OF_RATIO + 1)
"""

# MATERIAL PARAMETERS
THERMAL_K       = 385       # [W/mK]    - thermal conductivity of the material
FRICTION        = 0.012     # [-]       - friction factor

# GEOMETRICAL PARAMETERS
N_CHANNELS      = 20        # [-]       - number of channels in the combustion chamber
INTER_CHANNEL_T = 0.005     # [m]       - thickness of the wall separating the channels
LENGHT_CC       = 1.0       # [m]       - length of the combustion chamber
HEIHT_CHANNEL   = 0.003     # [m]       - height of the cooling channel
DI              = 0.70      # [m]       - inner diameter of the combustion chamber
T               = 0.01      # [m]       - thickness of the wall separating the combustion chamber from the coolant circuit

# HOT GAS PARAMETERS - INSIDE THE COMBUSTION CHAMBER (CC)
HOT_GAS_T_CC    = 2000   # [K]                   - temperature of the hot gases
HOT_GAS_CP_CC   = 6950      # [J/kgK]     REVIEW    - specific heat of the hot gases
HOT_GAS_MU_CC   = 1.17e-4   # [Pa.s]                - viscosity of the hot gases
HOT_GAS_K_CC    = 1.47566   # [W/mK]      REVIEW    - thermal conductivity of the hot gases
HOT_GAS_RHO_CC  = 9.1444    # [kg/m3]               - density of the hot gases
HOT_GAMMA_CC    = 1.1318    # [-]                   - specific heat ratio of the hot gases
HOT_SON_V_CC    = 1243.8    # [m/s]                 - speed of sound of the hot gases

# INITIAL ASSUMPTIONS
T_WI_INIT_ASSUMPTION = HOT_GAS_T_CC - 1000 # [K] - initial temperature assumption for inner combustion chamber wall
T_WO_INIT_ASSUMPTION = HOT_GAS_T_CC - 1200 # [K] - initial temperature assumption for outer combustion chamber wall

# COOLANT PROPERTIES
CONF_COOLANT_FILES = "conf/coolant"

# EUHAUST = 6.1 bar
