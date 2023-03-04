# Tools
import pyromat as pm
from math import pi


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
            section.get_exit_conditions()

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


    def get_ohmic_thermal_equivalences ():

        # NOTE: This is the thermal resistance of the coolant circuit
        R_cc = self.t / (self.k * pi * self.Di) # [K/W] - thermal resistance of the coolant circuit



    def get_exit_condition (self):

        

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
        e = self.t2 * self.h * (self.n-1)                           # [m2] - area lost to channel separation walls
        A = pi/4 * (pow(self.Di + self.h, 2) - pow(self.Di, 2))     # [m2] - cross-sectional area of the coolant circuit

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

        



# Main
def main ():
    combustion_chamber = CombustionChamber()
    combustion_chamber.execute()


if __name__ == "__main__":
    main()
