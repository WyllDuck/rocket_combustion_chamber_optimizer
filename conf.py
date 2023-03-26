""" GLOBAL PARAMETERS """

# COOLANT PARAMETERS
INLET_T         = 200       # [K]       - inlet temperature coolant
INLET_P         = 145       # [bar]     - inlet pressure coolant
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
MELTING_T       = 1357      # [K]       - melting temperature of the material - cooper

# GEOMETRICAL PARAMETERS
N_CHANNELS          = 200       # [-]       - number of channels in the combustion chamber
INTER_CHANNEL_T     = 0.001     # [m]       - thickness of the wall separating the channels
HEIHT_CHANNEL       = 0.01      # [m]       - height of the cooling channel
DI_CC               = 0.434405  # [m]       - inner diameter of the combustion chamber injection plate
DI_TH               = 0.2410    # [m]       - inner diameter of the combustion chamber throat
T                   = 0.005     # [m]       - thickness of the wall separating the combustion chamber from the coolant circuit

# INITIAL ASSUMPTIONS
T_WI_INIT_ASSUMPTION = 1300 # [K] - initial temperature assumption for inner combustion chamber wall
T_WO_INIT_ASSUMPTION = 1200 # [K] - initial temperature assumption for outer combustion chamber wall

# COOLANT PROPERTIES
CONF_COOLANT_FILES = "conf/coolant"
CONF_GAS_FILES     = "conf/gas"
GEOMETRY_FILE      = "geo/geometry.csv"

#CEA_FILES           = ["CEA/CH4 gaseous OF3_6 transport properties.txt", "C:/Users/marti/Desktop/rocket propulsion/rocket_combustion_chamber_optimizer/CEA/CH4 gaseous OF3_6 transport properties close to throat.txt"]
CEA_FILES           = ["CEA/CH4 gaseous OFvary transport properties.txt", "CEA/CH4_2.0_expansion_big.txt"]

# EUHAUST = 6.1 bar
