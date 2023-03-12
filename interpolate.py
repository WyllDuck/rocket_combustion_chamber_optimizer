# Tools
from scipy.interpolate import LinearNDInterpolator, interp1d
import pandas as pd
from math import pi

# Global variables
from conf import *


def get_interpolation_function_coolant_properties ():

    # function_cp  --> [kJ/kgK]  - function to get the specific heat capacity of the coolant
    # function_k   --> [mW/mK]   - function to get the thermal conductivity of the coolant
    # function_mu  --> [uPa*s]   - function to get the viscosity of the coolant
    # function_rho --> [kg/m3]   - function to get the density of the coolant

    # Open .csv file using pandas convert to numpy array
    data_coolant_cp     = pd.read_csv("{}_cp.csv".format(CONF_COOLANT_FILES), header=0).to_numpy()   # T, cp
    data_coolant_mu     = pd.read_csv("{}_mu.csv".format(CONF_COOLANT_FILES), header=0).to_numpy()   # T, p, mu
    data_coolant_k      = pd.read_csv("{}_k.csv".format(CONF_COOLANT_FILES), header=0).to_numpy()    # T, p, k
    data_coolant_rho    = pd.read_csv("{}_rho.csv".format(CONF_COOLANT_FILES), header=0).to_numpy()  # T, p, rho

    function_cp  = interp1d(data_coolant_cp[:, 0], data_coolant_cp[:, 1], fill_value=6.95)
    function_mu  = LinearNDInterpolator(data_coolant_mu[:, 0:2], data_coolant_mu[:, 2], fill_value=15)
    function_k   = LinearNDInterpolator(data_coolant_k[:, 0:2], data_coolant_k[:, 2], fill_value=50)
    function_rho = LinearNDInterpolator(data_coolant_rho[:, 0:2], data_coolant_rho[:, 2], fill_value=9.24)

    return function_cp, function_mu, function_k, function_rho


def get_interpolation_function_gas_properties ():

    # function_T        --> [K]         - function to get the temperature of the gas in the nozzle
    # function_cp       --> [kJ/kgK]    - function to get the specific heat capacity of the gas in the nozzle
    # function_k        --> [mW/mK]     - function to get the thermal conductivity of the gas in the nozzle
    # function_mu       --> [uPa*s]     - function to get the viscosity of the gas in the nozzle
    # function_rho      --> [kg/m3]     - function to get the density of the gas in the nozzle
    # function_gamma    --> [-]         - function to get the specific heat ratio of the gas in the nozzle
    # function_c        --> [m/s]       - function to get the speed of sound of the gas in the nozzle

    # Open .csv file using pandas convert to numpy array
    data_gas_T      = pd.read_csv("{}_T.csv".format(CONF_GAS_FILES), header=0).to_numpy()       # expansion_ratio
    data_gas_cp     = pd.read_csv("{}_cp.csv".format(CONF_GAS_FILES), header=0).to_numpy()      # expansion_ratio
    data_gas_mu     = pd.read_csv("{}_mu.csv".format(CONF_GAS_FILES), header=0).to_numpy()      # expansion_ratio
    data_gas_k      = pd.read_csv("{}_k.csv".format(CONF_GAS_FILES), header=0).to_numpy()       # expansion_ratio
    data_gas_rho    = pd.read_csv("{}_rho.csv".format(CONF_GAS_FILES), header=0).to_numpy()     # expansion_ratio
    data_gas_gamma  = pd.read_csv("{}_gamma.csv".format(CONF_GAS_FILES), header=0).to_numpy()   # expansion_ratio
    data_gas_c      = pd.read_csv("{}_c.csv".format(CONF_GAS_FILES), header=0).to_numpy()       # expansion_ratio

    # Interolation for the noozle region
    function_T_nozzle      = interp1d(data_gas_T[1:, 0], data_gas_T[1:, 1])
    function_cp_nozzle     = interp1d(data_gas_cp[1:, 0], data_gas_cp[1:, 1])
    function_mu_nozzle     = interp1d(data_gas_mu[1:, 0], data_gas_mu[1:, 1])
    function_k_nozzle      = interp1d(data_gas_k[1:, 0], data_gas_k[1:, 1])
    function_rho_nozzle    = interp1d(data_gas_rho[1:, 0], data_gas_rho[1:, 1])
    function_gamma_nozzle  = interp1d(data_gas_gamma[1:, 0], data_gas_gamma[1:, 1])
    function_c_nozzle      = interp1d(data_gas_c[1:, 0], data_gas_c[1:, 1])

    functions_nozzle =  [function_T_nozzle, function_cp_nozzle, function_mu_nozzle, function_k_nozzle, function_rho_nozzle, function_gamma_nozzle, function_c_nozzle]

    # Interpolation for the combustion chamber region
    expansion_ratio_cc = pow(DI_CC, 2) / pow(DI_TH, 2)
    
    # Change the expansion ratio to the combustion chamber value
    data_gas_T[0,0]         = expansion_ratio_cc
    data_gas_cp[0,0]        = expansion_ratio_cc
    data_gas_mu[0,0]        = expansion_ratio_cc
    data_gas_k[0,0]         = expansion_ratio_cc
    data_gas_rho[0,0]       = expansion_ratio_cc
    data_gas_gamma[0,0]     = expansion_ratio_cc
    data_gas_c[0,0]         = expansion_ratio_cc

    function_T_cc           = interp1d(data_gas_T[:2, 0], data_gas_T[:2, 1])
    function_cp_cc          = interp1d(data_gas_cp[:2, 0], data_gas_cp[:2, 1])
    function_mu_cc          = interp1d(data_gas_mu[:2, 0], data_gas_mu[:2, 1])
    function_k_cc           = interp1d(data_gas_k[:2, 0], data_gas_k[:2, 1])
    function_rho_cc         = interp1d(data_gas_rho[:2, 0], data_gas_rho[:2, 1])
    function_gamma_cc       = interp1d(data_gas_gamma[:2, 0], data_gas_gamma[:2, 1])
    function_c_cc           = interp1d(data_gas_c[:2, 0], data_gas_c[:2, 1])

    functions_cc =  [function_T_cc, function_cp_cc, function_mu_cc, function_k_cc, function_rho_cc, function_gamma_cc, function_c_cc]

    return functions_cc, functions_nozzle