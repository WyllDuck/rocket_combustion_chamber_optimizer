# Tools
import pyromat as pm
from math import pi, log


# FIXED GLOBAL PARAMETERS
INLET_T     = 200   # [K]       - inlet temperature coolant
INLET_p     = 130   # [bar]     - inlet pressure coolant

THERMAL_K   = 385   # [W/mK]    - thermal conductivity of the material
N_SECTIONS  = 20    # [º]       - number of cross-sections of the combustion chamber - simulation resolution


class CombustionChamber (object):

    def __init__ (self, dimensions, thermal_parameters_hot_gases, flow) -> None:
        
        # Dimension parameters
        Di, L, t, t2, n = dimensions # [m] - inner diameter, length, wall thickness, wall thickness (channels separating wall), number of channels
        
        self.Di = Di    # [m] - inner diameter
        self.L  = L     # [m] - length cc

        self.t  = t     # [m] - wall thickness (of the combustion chamber to coolant circuit)
        self.n  = n     # [º] - number of channels


        """ Section parameters """
        mdot, = flow # [kg/s] - mass flow rate

        Section.dx = self.L / N_SECTIONS    # [m]       - length of each cross-section (all sections have the same length)
        Section.mdot = mdot                 # [kg/s]    - mass flow rate (all sections have the same mass flow rate)
        Section.n = n                       # [º]       - number of channels (all sections have the same number of channels)
        Section.t2 = t2                     # [m]       - wall thickness (channels separating wall)


        """ Instance cross-sections of the combustion chamber """
        self.section = list(range(N_SECTIONS))
        self.section[0] = Section() 

        for i in range (1, N_SECTIONS):
            self.section[i] = Section(self.section[i-1])


        """ Thermal parameters - hot gas """
        T, cp, mu, k, rho, v = thermal_parameters_hot_gases # [J/kgK], [Pa s], [W/mK], [kg/m3], [m/s]

        self.T_gas      = T     # [K]       - temperature of the hot gas
        self.cp_gas     = cp    # [J/kgK]   - specific heat capacity of the gas
        self.mu_gas     = mu    # [Pa s]    - dynamic viscosity of the gas
        self.k_gas      = k     # [W/mK]    - thermal conductivity of the gas
        self.rho_gas    = rho   # [kg/m3]   - density of the gas
        self.v_gas      = v     # [m/s]     - velocity of the gas
      
    
    # execute the simulation
    def execute (self):

        for section in self.sections:
            section.get_exit_state()

        return 


    def get_convection_coefficient (self, T_wi):
        """ 
        Returns the convection coefficient of the hot-gas inside the combustion chamber
        --------------------------------
        in: 
            T_wi: [K] - temperature of the inner wall of the combustion chamber
        out:
            h_gas: [W/m2K] - convection coefficient of the hot-gas inside the combustion chamber
        """
        
        Re = self.rho_gas * self.v_gas * self.Di / self.mu_gas
        Pr = self.cp_gas * self.mu_gas / self.k_gas
        Nu = 0.0162 * pow(Re * Pr, 0.82) * pow(self.T_gas / T_wi, 0.57)

        return self.k_gas * Nu / self.Di # hot-gas convection coefficient



# Cross section of the combustion chamber 
class Section (object):

    # NOTE: These parameters are shared by all sections
    dx      = None  # [m]       - length of each cross-section
    mdot    = None  # [kg/s]    - mass flow rate - per section (not channel!)
    n       = None  # [º]       - number of channels
    t2      = None  # [m]       - wall thickness (channels separating wall)


    def __init__ (self, previous_section = None) -> None:
        
        # find gas in pyromat
        self.gas = pm.get('ig.CH4')

        self.pre_section = previous_section     # section before this one.
        
        # NOTE: If a previous section is provided then override the values, if not use CC_INLET conditions
        if not self.pre_section:        
            self.T_in   = self.pre_section.T_out
            self.p_in   = self.pre_section.p_out

        # INLET conditions
        else:
            self.T_in   = INLET_T # global variable
            self.p_in   = INLET_p # global variable

        """ placeholder calculations latter """
        # OUTLET conditions - placeholder calculation latter
        self.T_out = 0
        self.p_out = 0

        # Dimension parameters
        # NOTE: Height from the outer wall of the combustion chamber to the inner wall of the coolant circuit. 
        self.h  = 0     # [m2]      - average height per channel, height inlet to outlet.

        # Mass flow rate
        self.v    = 0   # [m/s]     - velocity

        # Thermal parameters - default values
        self.update_thermal_properties(self.T_in, self.p_in)

        # Filled in execute() / get_exit_state()
        self.R_co = None
        self.R_cc = None
        self.R_w   = None
        self.q    = None


    
    # Get the exit conditions of the section
    def execute (self):
        
        self.get_exit_state()
        
        
        pass

    # Function execute the analysis of the section - find outlet state of the section
    def get_exit_state (self, T_hot_gas): # [K] - temperature of the hot gas entering the combustion chamber

        """ PRESSURE ANALYSIS """
        self.p_out = self.p_in - self.get_pressure_loss()
        p_ = 0.5 * (self.p_out + self.p_in)


        """ TEMPERATURE ANALYSIS """
        T_out_assumption = self.T_in + 10 # [K] - initial guess for the outlet temperature

        # Iterate until the exit temperature is found
        while abs(T_out_assumption - T_out_cp) > 1e-3:
            
            # Determine the exit temperature
            R_co, R_cc, R_w = self.get_ohmic_thermal_equivalences()

            # Average temperature and pressure in the section
            T_ = 0.5 * (T_out_assumption + self.T_in)

            # Heat transfer rate analysis - ohmic resistance + cp method
            q           = (T_hot_gas - T_) / (R_co + R_cc + R_w)     # [W/m] - specific heat transfer rate
            Q           = q * self.dx * (pi * (self.Di + 2*self.t)) # [W]   - total heat transfer rate
            T_out_cp    = self.T_in + self.mdot / (self.cp * Q)     # [K]   - outlet temperature via the cp method

            # Get the new thermal parameters
            self.update_thermal_properties(T_, p_)

        # Save the ohmic resistance and the specific heat transfer rate
        self.R_co = R_co
        self.R_cc = R_cc
        self.R_w  = R_w
        self.q    = q

        return self.T_out, self.p_out # [K], [Pa]


    # Get ohmic equivalent thermal resistances
    def get_ohmic_thermal_equivalences (self):

        R_co = 1 / (self.h_co * (self.Di + self.t) * pi)               # [K/W] - thermal resistance of the coolant circuit
        R_cc = 1 / (self.h_cc * (self.Di) * pi)                        # [K/W] - thermal resistance of the combustion chamber
        R_w   = log((self.Di + self.t) / self.Di) / (2*pi*THERMAL_K)    # [K/W] - thermal resistance of the wall

        return R_co, R_cc, R_w

        
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
        e = self.t2 * self.h * (self.n-1)                                           # [m2] - area lost to channel separation walls
        A = pi/4 * (pow(self.Di + self.h + self.t, 2) - pow(self.Di + self.t, 2))   # [m2] - cross-sectional area of the coolant circuit

        return self.mdot / ((A-e) * self.rho) # [m/s] - average velocity of the coolant in the channels for this section


    def get_convection_coefficient (self, T_wo):
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

    
    def update_thermal_properties (self, T, p):

        self.cp     = pm.gas.Cp(self.gas, T=T, p=p)
        self.mu     = pm.gas.mu(self.gas, T=T, p=p)
        self.k      = pm.gas.k(self.gas, T=T, p=p)
        self.rho    = pm.gas.rho(self.gas, T=T, p=p)


    def get_pressure_loss (self):
        """ 
        Returns the pressure loss in the coolant circuit for this section
        --------------------------------
        in: 
            None
        out:
            dp: [Pa] - pressure loss in the coolant circuit for this section
        """
        
        return 0


    # returns keys temperature of the sections (inner, outer walls CC and the T_in, T_out)
    # Used to asses the temperature distribution in the combustion chamber, and the validity of the design
    # melting temperature of the wall.
    def get_section_temperatures (self):


        if not self.R_cc or not self.R_co or not self.R_w:
            print("ERROR: get_section_temperatures() - section not executed yet.")
            return -1, -1, -1, -1

        # Determine the temperature of the inner and outer walls of the combustion chamber
        self.T_wi = T_hot_gas - self.R_cc * self.q
        self.T_wo = self.T_wi - self.R_w * self.q

        return T_wi, T_wo, self.T_in, self.T_out




# Main
def main ():
    combustion_chamber = CombustionChamber()
    combustion_chamber.execute()


if __name__ == "__main__":
    main()
