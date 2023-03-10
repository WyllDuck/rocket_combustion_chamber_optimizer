# Tools
from math import log10, sqrt, pi

# Global variables
from conf import *


""" COMBUSTION CHAMBER CLASS """
class CombustionChamber (object):

    def __init__ (self, Di_nsection) -> None:
        
        self.Di_nsection = Di_nsection # [m] - inner diameter of the combustion chamber for each section

        # Hot-gas properties
        self.mdot       = MDOT              # [kg/s]    - mass flow rate of the hot-gas inside the combustion chamber
        
        self.T_gas      = HOT_GAS_T_CC      # [K]       - temperature of the hot gas
        self.cp_gas     = HOT_GAS_CP_CC     # [J/kgK]   - specific heat capacity of the gas
        self.mu_gas     = HOT_GAS_MU_CC     # [Pa s]    - dynamic viscosity of the gas
        self.k_gas      = HOT_GAS_K_CC      # [W/mK]    - thermal conductivity of the gas
        self.rho_gas    = HOT_GAS_RHO_CC    # [kg/m3]   - density of the gas
        self.gamma_gas  = HOT_GAMMA_CC      # [-]       - specific heat ratio of the gas
        self.c_gas      = HOT_SON_V_CC      # [m/s]     - speed of sound of the gas
        

        # filled-in after execution
        self.v_gas      = None              # [m/s]     - velocity of the gas
      

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
