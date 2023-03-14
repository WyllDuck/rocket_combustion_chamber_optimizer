# Tools
import numpy as np

# Other files
from section import Section
from combustion_chamber import CombustionChamber
from interpolate import get_interpolation_function_coolant_properties, get_interpolation_function_gas_properties

# Global variables
from conf import *


""" SOLVER CLASS """
class Solver (object):

    def __init__(self) -> None:
        
        """ Combustion chamber parameters """

        functions_cc, functions_nozzle = get_interpolation_function_gas_properties()

        function_T, function_cp, function_mu, function_k, function_rho, function_gamma, function_c = functions_cc

        CombustionChamber.function_T_cc     = function_T   
        CombustionChamber.function_cp_cc    = function_cp
        CombustionChamber.function_mu_cc    = function_mu
        CombustionChamber.function_k_cc     = function_k
        CombustionChamber.function_rho_cc   = function_rho
        CombustionChamber.function_gamma_cc = function_gamma
        CombustionChamber.function_c_cc     = function_c

        function_T, function_cp, function_mu, function_k, function_rho, function_gamma, function_c = functions_nozzle

        CombustionChamber.function_T_nozzle     = function_T
        CombustionChamber.function_cp_nozzle    = function_cp
        CombustionChamber.function_mu_nozzle    = function_mu
        CombustionChamber.function_k_nozzle     = function_k
        CombustionChamber.function_rho_nozzle   = function_rho
        CombustionChamber.function_gamma_nozzle = function_gamma
        CombustionChamber.function_c_nozzle     = function_c

        """ Section parameters """
        Section.mdot    = MDOT_COOLANT              # [kg/s]    - mass flow rate (all sections have the same mass flow rate)
        Section.n       = N_CHANNELS                # [º]       - number of channels (all sections have the same number of channels)
        Section.t2      = INTER_CHANNEL_T           # [m]       - wall thickness (channels separating wall)
        Section.t       = T                         # [m]       - wall thickness (combustion chamber wall)
        Section.f       = FRICTION                  # [-]       - friction factor
        Section.h       = HEIHT_CHANNEL             # [m]       - height of the channel

        """ Coolant properties """
        function_cp, function_mu, function_k, function_rho = get_interpolation_function_coolant_properties()

        Section.function_cp     = function_cp
        Section.function_mu     = function_mu
        Section.function_k      = function_k
        Section.function_rho    = function_rho


        # --- Instantiate the sections
        # Open GEOMETRY file .csv and read the values
        data = np.loadtxt(GEOMETRY_FILE, delimiter=',')

        Di_section = []                             # [m] - diameter of the section
        self.sections = [None] * (len(data) - 1)    # [-] - list of sections

        for i in range (0, len(data) - 1):

            # Average diameter of the section
            Di_n = 0.5*(data[i,1]+data[i+1,1])
            Di_section.append(Di_n)

            # Instantiate the section
            self.sections[i] = Section(Di_n, abs(data[i+1,0]-data[i,0]), data[i,2])

        # --- Instantiate the combustion chamber
        self.cc         = CombustionChamber(Di_section)


    # Solve the problem
    def solve (self):
        for i in range (len(self.sections)):
            if not self.sections[i].validate_geometry():
                print("INVALID GEOMETRY")
                break
            self.solve_section(i)


    # Get the convection coefficients inside the combustion chamber and in the coolant circuit
    def update_all_convection_coefficients_section (self, T_wi_, T_wo_, n): # [º] - number of the section

        s = self.sections[n] # [º] - number of the section
        
        s.h_cc = self.cc.get_cc_convection_coefficient_nsection(T_wi_, s.Di)    # [W/m2K]   - convection coefficient of the hot-gas inside the combustion chamber
        s.h_co = s.get_co_convection_coefficient(T_wo_)                         # [W/m2K]   - convection coefficient of the coolant inside the combustion chamber

        return 0


    # Function execute the analysis of the section - find outlet state of the section
    def solve_section (self, n): # [º] - number of the section

        s = self.sections[n]                # [º] - number of the section
        region = self.sections[n].region    # [-] - region of the section
        i = 0                               # [-] - iteration counter

        if n > 0: 
            s.set_inlet_state(self.sections[n-1].T_out, self.sections[n-1].p_out, self.sections[n-1])
        else:
            s.set_inlet_state(INLET_T, INLET_p)

        s.Dh = s.get_hydraulic_diameter()   # [m] - hydraulic diameter of the section

        # Initial assumptions
        s.T_out   = s.T_in + 10             # [K]   - initial guess for the outlet temperature
        s.T_wi_   = T_WI_INIT_ASSUMPTION    # [K]   - initial guess for the inner wall temperature
        s.T_wo_   = T_WO_INIT_ASSUMPTION    # [K]   - initial guess for the outer wall temperature
        s.p_out   = s.p_in                  # [Pa]  - initial guess for the pressure

        T_out_assumption    = 0
        T_wi_assumption     = 0
        T_wo_assumption     = 0
        p_out_assumption    = 0

        Section.print_class_values()
        self.cc.print_instance_values()

        # booting the iteration with the initial assumptions
        p_ = 0.5 * (s.p_out + s.p_in)           # average pressure in the section
        T_ = 0.5 * (s.T_out + s.T_in)           # average temperature in the section            
        s.update_thermal_properties(T_, p_)     # [K], [Pa] - temperature, pressure - update the thermal properties of the coolant

        # Set the cc thermal properties for the section
        exp_ratio = pow(s.Di, 2) / pow(DI_TH, 2)
        self.cc.get_thermal_properties_gas(region, exp_ratio)
        self.cc.v_gas   = self.cc.get_velocity_hot_gases()      # [m/s] - velocity of the hot-gas inside the combustion chamber

        # Iterate until the exit temperature is found
        while sum(self.check_assumptions_section(n, (T_out_assumption, T_wi_assumption, T_wo_assumption, p_out_assumption))) > 1e-3:

            # Update the temperature assumptions - saved to latter check if the assumptions are correct
            T_out_assumption    = s.T_out
            T_wi_assumption     = s.T_wi_
            T_wo_assumption     = s.T_wo_
            p_out_assumption    = s.p_out

            s.p_out = s.p_in - s.get_pressure_loss()*0.00001    # [bar] - pressure loss in the section

            # Update thermal properties based on the new assumptions
            p_ = 0.5 * (s.p_out + s.p_in)           # average pressure in the section
            T_ = 0.5 * (s.T_out + s.T_in)           # average temperature in the section            
            s.update_thermal_properties(T_, p_)     # [K], [Pa] - temperature, pressure - update the thermal properties of the coolant

            s.v = s.get_velocity_coolant()          # [m/s] - velocity of the coolant in the section
            self.update_all_convection_coefficients_section(s.T_wi_, s.T_wo_, n)
            R_co, R_cc, R_w = s.get_ohmic_thermal_equivalences()

            # Heat transfer rate analysis - ohmic resistance + cp method
            Q           = (self.cc.T_gas - T_) / (R_co + R_cc + R_w)    # [W]   - heat transfer rate
            s.T_out     = s.T_in + Q / (s.cp * s.mdot)                  # [K]   - outlet temperature via the cp method

            # Save the ohmic resistance and the heat transfer rate
            s.R_co = R_co
            s.R_cc = R_cc
            s.R_w  = R_w
            s.Q    = Q

            print("--------------------------------------------")
            s.print_instance_values()

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
