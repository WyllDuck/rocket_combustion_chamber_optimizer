# Tools
from math import pi, log, sqrt, tanh

# Global variables
from conf import *


""" SECTION CLASS """
class Section (object):

    # NOTE: These parameters are shared by all sections
    mdot    = None  # [kg/s]    - mass flow rate - per section (not channel!)
    n       = None  # [º]       - number of channels
    t2      = None  # [m]       - wall thickness (channels separating wall)
    t       = None  # [m]       - wall thickness (combustion chamber wall)
    f       = None  # [-]       - friction factor
    h       = None  # [m]       - height of the cooling channel

    # coolant interpolation functions
    function_cp     = None
    function_mu     = None
    function_k      = None
    function_rho    = None


    def __init__ (self, Di, dx, region) -> None:

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

        self.T_wi_  = 0 # [K]       - temperature of the inner wall of the combustion chamber
        self.T_wo_  = 0 # [K]       - temperature of the outer wall of the combustion chamber

        self.Dh     = 0         # [m]   - hydraulic diameter of the cooling channels
        self.Di     = Di        # [m]   - inner diameter of the combustion chamber
        self.dx     = dx        # [m]   - length of the cross-section
        self.region = region    # [-]   - region of the rocket engine, inside the CC or in the nozzle

        self.thickness = 0      # [m]   - thickness of the wall


    # Get the hydraulic diameter of the coolant circuit
    def get_hydraulic_diameter (self):
        
        # Area of the coolant circuit cross-section
        e = self.t2 * self.h * self.n                                                       # [m2] - area lost to channel separation walls
        At = pi/4 * (pow(self.Di + 2*(self.h + self.t), 2) - pow(self.Di + 2*self.t, 2))    # [m2] - cross-sectional area of the coolant circuit
        A = (At - e)/self.n                                                                 # [m2] - subtract the area lost to channel separation walls + divide by number of cooling channels

        Pi = (pi*(self.Di + 2*self.t) - self.n*self.t2)/self.n                              # [m]  - length edge inner perimeter of the coolant channel
        Po = (pi*(self.Di + 2*(self.h + self.t)) - self.n*self.t2)/self.n                   # [m]  - length edge outer perimeter of the coolant channel
        P = Pi + Po + 2*self.h                                                              # [m]  - total perimeter of the coolant channel

        print("AR_channel: \t", HEIGHT_CHANNEL/(0.5*(Pi+Po)))

        return 4*A/P


    # Set Inlet State
    def set_inlet_state (self, T_in, p_in, previous_section = None):

        # Set initial assumptions for wall temperatures - speeds up the solver
        if previous_section is not None:
            self.T_wi_  = previous_section.T_wi_
            self.T_wo_  = previous_section.T_wo_

        self.T_in   = T_in
        self.p_in   = p_in

        # boot start solver
        self.p_out  = self.p_in   # [Pa]  - outlet pressure
        self.T_out  = self.T_in   # [K]   - outlet temperature


    # Check if the geometry is valid
    def validate_geometry (self):
        return self.n * self.t2 < (self.Di + 2*self.t) * pi


    # print all instance values
    def print_instance_values (self):
        
        print("--- coolant_properties ---")
        print("cp: \t", self.cp)
        print("mu: \t", self.mu)
        print("k: \t", self.k)
        print("rho: \t", self.rho)
        print()

        print("--- thermal_resistances ---")
        print("R_co: \t", self.R_co)
        print("R_cc: \t", self.R_cc)
        print("R_w: \t", self.R_w)
        print()

        print("--- convection_coefficients ---")
        print("h_cc: \t", self.h_cc)
        print("h_co: \t", self.h_co)
        print()

        print("--- other ---")
        print("v: \t", self.v)
        print("p_out: \t", self.p_out)
        print("Q: \t", self.Q)
        print("Di: \t", self.Di)
        print("Dh: \t", self.Dh)
        print("region: \t", self.region)
        print("thickness: \t", self.thickness)
        print()

        print("--- temperatures ---")
        print("T_out: \t", self.T_out)
        print("T_in: \t", self.T_out)
        print("T_wi_: \t", self.T_wi_)
        print("T_wo_: \t", self.T_wo_)


    # print all class values
    @staticmethod
    def print_class_values():

        print("--- class values ---")
        print("mdot: \t", Section.mdot)
        print("n: \t", Section.n)
        print("t2: \t", Section.t2)
        print("t: \t", Section.t)
        print("f: \t", Section.f)
        print("h: \t", Section.h)
        print()

    def get_areas_and_perimeters_heat_transfer (self):
            
        Pwi = self.Di*pi                                            # [m]  - inner perimeter combustion chamber
        Pwo = (self.Di + 2*self.t) * pi - self.n * self.t2          # [m]  - interface wall perimeter from the coolant circuit with combustion chamber

        Awi = Pwi * self.dx                                         # [m2] - inner area of the combustion chamber
        Awo = Pwo * self.dx                                         # [m2] - interface area of the wall from the coolant circuit with combustion chamber

        Af  = self.h * self.dx * 2                                  # [m2] - Heat tranfer area for 1 fin

        # Values for the fin modeling - obtain fin parameter 'm'
        Pf  = self.dx * 2                                           # [m]  - Heat transfer perimeter
        Ac  = self.dx * self.t2                                     # [m2] - Cross sectional area fin

        return (Awi, Awo, Af), (Pwi, Pwo), (Pf, Ac)


    # Get ohmic equivalent thermal resistances
    def get_ohmic_thermal_equivalences (self):

        (Awi, Awo, Af), (Pwi, Pwo), (Pf, Ac) = self.get_areas_and_perimeters_heat_transfer()
        
        R_co = 1 / (self.h_co * Awo)                                                    # [K/W] - thermal resistance of the coolant circuit
        R_cc = 1 / (self.h_cc * Awi)                                                    # [K/W] - thermal resistance of the combustion chamber
        R_w  = log((self.Di + 2*self.t) / self.Di) / (2*pi*THERMAL_K*self.dx)           # [K/W] - thermal resistance of the wall

        # Fin Modeling
        m   = sqrt((self.h_co * Pf)/(THERMAL_K * Ac))
        eff = tanh(m*self.h) / (m*self.h)
        R_f = 1 / (eff * Af * self.h_co)

        print("eff: ", eff)
        print("m: ", m)

        print("R_co: ", R_co)
        R_co = pow((self.n/R_f + 1/R_co), -1)
        print("R_co: ", R_co)

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
        e = self.t2 * self.h * self.n                                                   # [m2] - area lost to channel separation walls
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
        
        Re = self.rho * self.v * self.Dh / self.mu
        Pr = self.cp * self.mu / self.k
        Nu = 0.023 * pow(Re, 0.8) * pow(Pr, 0.4) * pow(0.5*(self.T_in+self.T_out) / T_wo, 0.57)

        return self.k * Nu / self.Dh # hot-gas convection coefficient


    # Update the thermal properties of the coolant
    def update_thermal_properties (self, _T, _p):

        self.cp     = self.function_cp(_T)       * 1e3           # [J/kgK] - specific heat capacity of the coolant
        self.mu     = self.function_mu(_T, _p)   * 1e-6          # [Pa s] - dynamic viscosity of the coolant
        self.k      = self.function_k(_T, _p)    * 1e-3          # [W/mK] - thermal conductivity of the coolant
        self.rho    = self.function_rho(_T, _p)                  # [kg/m3] - density of the coolant


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
        
        return self.f * (self.rho * pow(self.v, 2) / 2) * self.dx / self.Dh  # [Pa] - pressure loss in the coolant circuit for this section
    

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