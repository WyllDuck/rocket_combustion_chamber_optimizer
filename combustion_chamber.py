# Tools
from math import log10, sqrt, pi

# Global variables
from conf import *


""" COMBUSTION CHAMBER CLASS """
class CombustionChamber (object):

    # Gas Thermal Properties
    # from injection plate to combustion chamber throat
    function_T_cc            = None      # [K]       - temperature of the hot gas inside the combustion chamber
    function_cp_cc           = None      # [J/kgK]   - specific heat capacity of the gas inside the combustion chamber
    function_mu_cc           = None      # [Pa s]    - dynamic viscosity of the gas inside the combustion chamber
    function_k_cc            = None      # [W/mK]    - thermal conductivity of the gas inside the combustion chamber
    function_rho_cc          = None      # [kg/m3]   - density of the gas inside the combustion chamber
    function_gamma_cc        = None      # [-]       - specific heat ratio of the gas inside the combustion chamber
    function_c_cc            = None      # [m/s]     - speed of sound of the gas inside the combustion chamber
    
    # from combustion chamber throat to nozzle exit
    function_T_nozzle        = None      # [K]       - temperature of the hot gas at the nozzle positions of the combustion chamber
    function_cp_nozzle       = None      # [J/kgK]   - specific heat capacity of the gas at the nozzle positions of the combustion chamber
    function_mu_nozzle       = None      # [Pa s]    - dynamic viscosity of the gas at the nozzle positions of the combustion chamber
    function_k_nozzle        = None      # [W/mK]    - thermal conductivity of the gas at the nozzle positions of the combustion chamber
    function_rho_nozzle      = None      # [kg/m3]   - density of the gas at the nozzle positions of the combustion chamber
    function_gamma_nozzle    = None      # [-]       - specific heat ratio of the gas at the nozzle positions of the combustion chamber
    function_c_nozzle        = None      # [m/s]     - speed of sound of the gas at the nozzle positions of the combustion chamber


    def __init__ (self, Di_nsection) -> None:
        
        self.Di_nsection    = Di_nsection       # f(x) [m] - inner diameter of the combustion chamber for each section based on CC position
        self.mdot           = MDOT              # [kg/s]    - mass flow rate of the hot-gas inside the combustion chamber
    
        # filled-in after execution
        self.v_gas          = None              # [m/s]     - velocity of the gas

        # Gas Thermal Properties
        self.T_gas          =  None      # [K]       - temperature of the hot gas inside the combustion chamber
        self.cp_gas         =  None      # [J/kgK]   - specific heat capacity of the gas inside the combustion chamber
        self.mu_gas         =  None      # [Pa s]    - dynamic viscosity of the gas inside the combustion chamber
        self.k_gas          =  None      # [W/mK]    - thermal conductivity of the gas inside the combustion chamber
        self.rho_gas        =  None      # [kg/m3]   - density of the gas inside the combustion chamber
        self.gamma_gas      =  None      # [-]       - specific heat ratio of the gas inside the combustion chamber
        self.c_gas          =  None      # [m/s]     - speed of sound of the gas inside the combustion chamber
      

    def get_thermal_properties_gas (self, region, exp_ratio):
        """ 
        Returns the thermal properties of the hot-gas inside the combustion chamber
        --------------------------------
        in: 
            - region:       boolean value indicating the region of study: inside the combustion chamber or at the nozzle exit (1, 2)
            - exp_ratio:    expansion ratio
        out:
            -
        """
        
        # NOTE: UNITS NEED TO BE CONVERTED TO SI UNITS
        if region == 1:
            self.T_gas      = self.function_T_cc(exp_ratio)
            self.cp_gas     = self.function_cp_cc(exp_ratio)    * 1e3
            self.mu_gas     = self.function_mu_cc(exp_ratio)    * 1e-4
            self.k_gas      = self.function_k_cc(exp_ratio)     * 1e-1
            self.rho_gas    = self.function_rho_cc(exp_ratio)
            self.gamma_gas  = self.function_gamma_cc(exp_ratio)
            self.c_gas      = self.function_c_cc(exp_ratio)

        elif region == 2:
            self.T_gas      = self.function_T_nozzle(exp_ratio)
            self.cp_gas     = self.function_cp_nozzle(exp_ratio)    * 1e3
            self.mu_gas     = self.function_mu_nozzle(exp_ratio)    * 1e-4
            self.k_gas      = self.function_k_nozzle(exp_ratio)     * 1e-1
            self.rho_gas    = self.function_rho_nozzle(exp_ratio)
            self.gamma_gas  = self.function_gamma_nozzle(exp_ratio)
            self.c_gas      = self.function_c_nozzle(exp_ratio)
        

    def get_velocity_hot_gases (self):
        """ 
        Returns the velocity of the hot-gas inside the combustion chamber
        --------------------------------
        in: 
            -
        out:
            v_gas:  [m/s]   - velocity of the hot-gas inside the combustion chamber
        """
        
        # NOTE: The velocity is calculated for the last section because that is the injection plate section
        v_gas = self.mdot / (self.rho_gas * pi * pow(self.Di_nsection[-1], 2) / 4)
        return v_gas
    

    # print the values of the instance
    def print_instance_values (self):

        print("T_gas: \t\t{}".format(self.T_gas))
        print("cp_gas: \t{}".format(self.cp_gas))
        print("mu_gas: \t{}".format(self.mu_gas))
        print("k_gas: \t\t{}".format(self.k_gas))
        print("rho_gas: \t{}".format(self.rho_gas))
        print("v_gas: \t\t{}".format(self.v_gas))
        print("gamma_gas: \t{}".format(self.gamma_gas))
        print("c_gas: \t\t{}".format(self.c_gas))

        print()


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
        
        # OPTION 1
        if True:
            T_aw = self.T_gas * (1 + pow(Pr, 1/3) * (self.gamma_gas - 1) / 2 * pow(self.v_gas / self.c_gas, 2))
            Nu = 0.0162 * pow(Re * Pr, 0.82) * pow(T_aw / T_wi, 0.35)
            print("T_aw: \t{}".format(T_aw))

        # OPTION 2 
        else:
            epsilon = (1.82*log10(Re) - 1.64)**(-2)
            Nu = (epsilon/8) * (Re - 1000) * Pr / (1 + 12.7 * sqrt(epsilon/8) * (pow(Pr, 2/3) - 1)) * (1 + pow((Di / LENGHT_CC), 2/3))
        
        print("Re: \t{}".format(Re))
        print("Pr: \t{}".format(Pr))
        print("Nu: \t{}".format(Nu))

        return self.k_gas * Nu / Di # hot-gas convection coefficient
