# Tools
from math import log10, sqrt

# Global variables
from conf import *


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
      

    # print the values of the instance
    def print_instance_values (self):

        print("T_gas: \t\t{}".format(self.T_gas))
        print("cp_gas: \t{}".format(self.cp_gas))
        print("mu_gas: \t{}".format(self.mu_gas))
        print("k_gas: \t\t{}".format(self.k_gas))
        print("rho_gas: \t{}".format(self.rho_gas))
        print("v_gas: \t\t{}".format(self.v_gas))
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
        #Nu = 0.0162 * pow(Re * Pr, 0.82) * pow(self.T_gas / T_wi, 0.35)

        epsilon = (1.82*log10(Re) - 1.64)**(-2)
        Nu = (epsilon/8) * (Re - 1000) * Pr / (1 + 12.7 * sqrt(epsilon/8) * (pow(Pr, 2/3) - 1)) * (1 + pow((Di / LENGHT_CC), 2/3))

        print("Re: \t{}".format(Re))
        print("Pr: \t{}".format(Pr))
        print("Nu: \t{}".format(Nu))

        return self.k_gas * Nu / Di # hot-gas convection coefficient
