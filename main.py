# Tools
from math import pi, log
from scipy.interpolate import LinearNDInterpolator, interp1d
import pandas as pd


# FIXED GLOBAL PARAMETERS
N_SECTIONS      = 20    # [º]       - number of cross-sections of the combustion chamber - simulation resolution

# COOLANT PARAMETERS
INLET_T         = 200   # [K]       - inlet temperature coolant
INLET_p         = 130   # [bar]     - inlet pressure coolant
MDOT_COOLANT    = 0.1   # [kg/s]    - mass flow rate of the coolant

# MATERIAL PARAMETERS
THERMAL_K       = 385   # [W/mK]    - thermal conductivity of the material
FRICTION        = 0.01  # [-]       - friction factor

# GEOMETRICAL PARAMETERS
N_CHANNELS      = 10    # [-]       - number of channels in the combustion chamber
INTER_CHANNEL_T = 0.001 # [m]       - thickness of the wall separating the channels
LENGHT_CC       = 0.5   # [m]       - length of the combustion chamber
HEIHT_CHANNEL   = 0.1   # [m]       - height of the cooling channel
DI              = 0.1   # [m]       - inner diameter of the combustion chamber
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
T_WI_INIT_ASSUMPTION = HOT_GAS_T - 10 # [K] - initial temperature assumption for inner combustion chamber wall
T_WO_INIT_ASSUMPTION = HOT_GAS_T - 20 # [K] - initial temperature assumption for outer combustion chamber wall


""" SOLVER CLASS """
class Solver (object):

    def __init__(self) -> None:
        
        self.cc         = CombustionChamber(THERMAL_PARA_HOT_GASES, [DI] * N_SECTIONS)
        
        """ Section parameters """
        Section.dx      = LENGHT_CC / N_SECTIONS    # [m]       - length of each cross-section (all sections have the same length)
        Section.mdot    = MDOT_COOLANT              # [kg/s]    - mass flow rate (all sections have the same mass flow rate)
        Section.n       = N_CHANNELS                # [º]       - number of channels (all sections have the same number of channels)
        Section.t2      = INTER_CHANNEL_T           # [m]       - wall thickness (channels separating wall)
        Section.t       = T                         # [m]       - wall thickness (combustion chamber wall)
        Section.f       = FRICTION                  # [-]       - friction factor
        Section.h       = HEIHT_CHANNEL             # [m]       - height of the channel

        # Instantiate the sections
        self.sections = [None] * N_SECTIONS
        self.sections[0] = Section(self.cc.Di_nsection[0]) 

        for i in range (1, N_SECTIONS):
            self.sections[i] = Section(self.cc.Di_nsection[i], self.sections[i-1])


    # Solve the problem
    def solve (self):
        for i in range (N_SECTIONS):
            self.solve_section(i)


    # Get the convection coefficients inside the combustion chamber and in the coolant circuit
    def update_all_convection_coefficients_section (self, T_wi_, T_wo_, n): # [º] - number of the section

        s = self.sections[n] # [º] - number of the section
        
        s.h_cc = self.cc.get_cc_convection_coefficient_nsection(T_wi_, s.Di)    # [W/m2K]   - convection coefficient of the hot-gas inside the combustion chamber
        s.h_co = s.get_co_convection_coefficient(T_wo_)                         # [W/m2K]   - convection coefficient of the coolant inside the combustion chamber

        return 0


    # Function execute the analysis of the section - find outlet state of the section
    def solve_section (self, n): # [º] - number of the section

        s = self.sections[n] # [º] - number of the section
        i = 0                # [-] - iteration counter

        # Initial assumptions
        s.T_out   = s.T_in + 10             # [K]   - initial guess for the outlet temperature
        s.T_wi_   = T_WI_INIT_ASSUMPTION    # [K]   - initial guess for the inner wall temperature
        s.T_wo_   = T_WO_INIT_ASSUMPTION    # [K]   - initial guess for the outer wall temperature
        s.p_out   = s.p_in                  # [Pa]  - initial guess for the pressure

        T_out_assumption    = 0
        T_wi_assumption     = 0
        T_wo_assumption     = 0
        p_out_assumption    = 0
        

        # booting the iteration with the initial assumptions
        p_ = 0.5 * (s.p_out + s.p_in)           # average pressure in the section
        T_ = 0.5 * (s.T_out + s.T_in)           # average temperature in the section            
        s.update_thermal_properties(T_, p_)     # [K], [Pa] - temperature, pressure - update the thermal properties of the coolant


        # Iterate until the exit temperature is found
        while sum(self.check_assumptions_section(n, (T_out_assumption, T_wi_assumption, T_wo_assumption, p_out_assumption))) > 1e-3:

            # Update the temperature assumptions - saved to latter check if the assumptions are correct
            T_out_assumption    = s.T_out
            T_wi_assumption     = s.T_wi_
            T_wo_assumption     = s.T_wo_
            p_out_assumption    = s.p_out

            s.p_out = s.p_in - s.get_pressure_loss() # [Pa] - pressure loss in the section

            # Update thermal properties based on the new assumptions
            p_ = 0.5 * (s.p_out + s.p_in)           # average pressure in the section
            T_ = 0.5 * (s.T_out + s.T_in)           # average temperature in the section            
            s.update_thermal_properties(T_, p_)     # [K], [Pa] - temperature, pressure - update the thermal properties of the coolant

            self.update_all_convection_coefficients_section(s.T_wi_, s.T_wo_, n)
            R_co, R_cc, R_w = s.get_ohmic_thermal_equivalences()

            # Heat transfer rate analysis - ohmic resistance + cp method
            q           = (self.cc.T_gas - T_) / (R_co + R_cc + R_w)    # [W/m] - specific heat transfer rate
            Q           = q * s.dx * (pi * (s.Di + 2*s.t))              # [W]   - total heat transfer rate
            s.T_out     = s.T_in + s.mdot / (s.cp * Q)                  # [K]   - outlet temperature via the cp method

            # Save the ohmic resistance and the specific heat transfer rate
            s.R_co = R_co
            s.R_cc = R_cc
            s.R_w  = R_w
            s.Q    = Q

            # Check other assumptions (wall temperature update)
            s.get_section_temperatures(self.cc.T_gas)           # [K]   - get the wall temperatures - inner and outer

            i += 1  # update the iteration counter
            print("section {} - iteration {} - T_out = {} K".format(n, i, s.T_out))

        return s.T_out, s.p_out # [K], [Pa]
    

    # check the temperature assumptions validity
    def check_assumptions_section (self, n, assumptions): # [º] - number of the section

        s = self.sections[n]
        T_out_assumption, T_wi_assumption, T_wo_assumption, p_out_assumption = assumptions

        return abs(s.T_wi_ - T_wi_assumption), abs(s.T_wo_ - T_wo_assumption), abs(s.T_out - T_out_assumption), abs(s.p_out - p_out_assumption)


""" COMBUSTION CHAMBER CLASS """
class CombustionChamber (object):

    def __init__ (self, thermal_parameters_hot_gases, Di_nsection) -> None:
        
        self.Di_nsection = Di_nsection # [m] - inner diameter of the combustion chamber for each section

        """ Thermal parameters - hot gas """
        T, cp, mu, k, rho, v = thermal_parameters_hot_gases # [J/kgK], [Pa s], [W/mK], [kg/m3], [m/s]

        self.T_gas      = T     # [K]       - temperature of the hot gas
        self.cp_gas     = cp    # [J/kgK]   - specific heat capacity of the gas
        self.mu_gas     = mu    # [Pa s]    - dynamic viscosity of the gas
        self.k_gas      = k     # [W/mK]    - thermal conductivity of the gas
        self.rho_gas    = rho   # [kg/m3]   - density of the gas
        self.v_gas      = v     # [m/s]     - velocity of the gas
      

    # Get the convection coefficient of the hot-gas inside the combustion chamber
    def get_cc_convection_coefficient_nsection (self, T_wi, Di): 
        """ 
        Returns the convection coefficient of the hot-gas inside the combustion chamber
        --------------------------------
        in: 
            T_wi:   [K] - temperature of the inner wall of the combustion chamber
            Di:     [m] - inner diameter of the combustion chamber
        out:
            h_gas:  [W/m2K] - convection coefficient of the hot-gas inside the combustion chamber
        """
        
        Re = self.rho_gas * self.v_gas * Di / self.mu_gas
        Pr = self.cp_gas * self.mu_gas / self.k_gas
        Nu = 0.0162 * pow(Re * Pr, 0.82) * pow(self.T_gas / T_wi, 0.57)

        return self.k_gas * Nu / Di # hot-gas convection coefficient


""" SECTION CLASS """
class Section (object):

    # NOTE: These parameters are shared by all sections
    dx      = None  # [m]       - length of each cross-section
    mdot    = None  # [kg/s]    - mass flow rate - per section (not channel!)
    n       = None  # [º]       - number of channels
    t2      = None  # [m]       - wall thickness (channels separating wall)
    t       = None  # [m]       - wall thickness (combustion chamber wall)
    f       = None  # [-]       - friction factor
    h       = None  # [m]       - height of the cooling channel

    gas = pm.get('ig.CH4') # [pyromat] - gas in the combustion chamber


    def __init__ (self, Di, previous_section = None) -> None:
        
        # find gas in pyromat
        self.pre_section = previous_section     # section before this one.
        
        # NOTE: If a previous section is provided then override the values, if not use CC_INLET conditions
        if self.pre_section:        
            self.T_in   = self.pre_section.T_out
            self.p_in   = self.pre_section.p_out

        # INLET conditions
        else:
            self.T_in   = INLET_T # global variable
            self.p_in   = INLET_p # global variable

        # Filled-in after execution

        self.cp     = 0 # [J/kgK]   - specific heat capacity of the coolant
        self.mu     = 0 # [Pa s]    - dynamic viscosity of the coolant
        self.k      = 0 # [W/mK]    - thermal conductivity of the coolant
        self.rho    = 0 # [kg/m3]   - density of the coolant

        self.R_co   = 0 # [K/W]     - thermal resistance of the coolant circuit convection
        self.R_cc   = 0 # [K/W]     - thermal resistance of the combustion chamber coonvection
        self.R_w    = 0 # [K/W]     - thermal resistance of the wall
        self.Q      = 0 # [W]       - heat transfer rate

        self.h_cc   = 0 # [W/m2K]   - convection coefficient of the hot-gas inside the combustion chamber
        self.h_co   = 0 # [W/m2K]   - convection coefficient of the coolant

        self.v      = 0 # [m/s]     - velocity

        self.p_out  = self.p_in   # [Pa]  - outlet pressure
        self.T_out  = self.T_in   # [K]   - outlet temperature

        self.T_wi_  = 0 # [K]       - temperature of the inner wall of the combustion chamber
        self.T_wo_  = 0 # [K]       - temperature of the outer wall of the combustion chamber

        self.Di     = Di # [m]       - inner diameter of the combustion chamber


    # Get ohmic equivalent thermal resistances
    def get_ohmic_thermal_equivalences (self):

        R_co = 1 / (self.h_co * (self.Di + 2*self.t) * pi)                # [K/W] - thermal resistance of the coolant circuit
        R_cc = 1 / (self.h_cc * (self.Di) * pi)                         # [K/W] - thermal resistance of the combustion chamber
        R_w  = log((self.Di + 2*self.t) / self.Di) / (2*pi*THERMAL_K)     # [K/W] - thermal resistance of the wall

        return R_co, R_cc, R_w


    # Get the velocity of the coolant in the channels        
    def get_velocity_coolant (self):

        """
        Returns the avergare velocity of the coolant in the channels for this section.
        --------------------------------
        in: 
            None
        out:
            v: [m/s] - average velocity of the coolant in the channels for this section
        """

        # Area of the coolant circuit cross-section
        e = self.t2 * self.h * (self.n-1)                                               # [m2] - area lost to channel separation walls
        A = pi/4 * (pow(self.Di + 2*(self.h + self.t), 2) - pow(self.Di + 2*self.t, 2)) # [m2] - cross-sectional area of the coolant circuit

        return self.mdot / ((A-e) * self.rho) # [m/s] - average velocity of the coolant in the channels for this section


    # Get the convection coefficient of the coolant inside the coolant circuit
    def get_co_convection_coefficient (self, T_wo):
        """ 
        Returns the convection coefficient of the hot-gas inside the combustion chamber
        --------------------------------
        in: 
            T_wo: [K] - temperature of the outer wall of the combustion chamber, coolant circuit
        out:
            h: [W/m2K] - convection coefficient of the coolant inside the coolant circuit
        """
        
        Re = self.rho * self.v * self.Di / self.mu
        Pr = self.cp * self.mu / self.k
        Nu = 0.023 * pow(Re, 0.8) * pow(Pr, 0.4) * pow(self.T_in / T_wo, 0.57)

        return self.k * Nu / self.Di # hot-gas convection coefficient


    # Update the thermal properties of the coolant
    def update_thermal_properties (self, _T, _p):

        self.cp     = FUNCTION_CP(_T)       * 1e3           # [J/kgK] - specific heat capacity of the coolant
        self.mu     = FUNCTION_MU(_T, _p)   * 1e-6          # [Pa s] - dynamic viscosity of the coolant
        self.k      = FUNCTION_K(_T, _p)    * 1e-3          # [W/mK] - thermal conductivity of the coolant
        self.rho    = FUNCTION_RHO(_T, _p)                  # [kg/m3] - density of the coolant


    # Get the pressure loss in the coolant circuit for this section
    def get_pressure_loss (self):
        """ 
        Returns the pressure loss in the coolant circuit for this section
        --------------------------------
        in: 
            None
        out:
            dp: [Pa] - pressure loss in the coolant circuit for this section
        """
        
        return self.f * (self.rho * pow(self.v, 2) / 2) * self.dx / self.h  # [Pa] - pressure loss in the coolant circuit for this section
    

    # Get the temperature of the inner and outer walls of the combustion chamber, and the inlet and outlet of the coolant circuit
    def get_section_temperatures (self, T_hot_gas):

        if not self.R_cc or not self.R_co or not self.R_w:
            print("ERROR: get_section_temperatures() - section not executed yet.")
            return -1, -1, -1, -1

        # Determine the temperature of the inner and outer walls of the combustion chamber
        # NOTE: The _ marks the variable is an average value over the whole section for this particular value.
        self.T_wi_ = T_hot_gas  - self.R_cc * self.Q
        self.T_wo_ = self.T_wi_ - self.R_w  * self.Q

        return self.T_wi_, self.T_wo_, self.T_in, self.T_out


def get_interpolation_function_for_mu_and_rho ():

    global FUNCTION_CP  # [kJ/kgK]  - function to get the specific heat capacity of the coolant
    global FUNCTION_K   # [mW/mK]   - function to get the thermal conductivity of the coolant
    global FUNCTION_MU  # [uPa*s]   - function to get the viscosity of the coolant
    global FUNCTION_RHO # [kg/m3]   - function to get the density of the coolant

    # Open .csv file using pandas convert to numpy array
    data_coolant_cp     = pd.read_csv("conf/coolant_cp.csv", header=0).to_numpy()   # T, cp
    data_coolant_mu     = pd.read_csv("conf/coolant_mu.csv", header=0).to_numpy()   # T, p, mu
    data_coolant_k      = pd.read_csv("conf/coolant_k.csv", header=0).to_numpy()    # T, p, k
    data_coolant_rho    = pd.read_csv("conf/coolant_rho.csv", header=0).to_numpy()  # T, p, rho

    FUNCTION_CP  = interp1d(data_coolant_cp[:, 0], data_coolant_cp[:, 1])
    FUNCTION_MU  = LinearNDInterpolator(data_coolant_mu[:, 0:2], data_coolant_mu[:, 2])
    FUNCTION_K   = LinearNDInterpolator(data_coolant_k[:, 0:2], data_coolant_k[:, 2])
    FUNCTION_RHO = LinearNDInterpolator(data_coolant_rho[:, 0:2], data_coolant_rho[:, 2])

    return 0


# Main
def main ():

    # Get the interpolation functions for the coolant
    get_interpolation_function_for_mu_and_rho()

    solver = Solver()
    solver.solve()


if __name__ == "__main__":
    main()
