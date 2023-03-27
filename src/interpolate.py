# Tools
from scipy.interpolate import LinearNDInterpolator, interp1d
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

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

    function_cp  = interp1d(data_coolant_cp[:, 0], data_coolant_cp[:, 1])
    function_mu  = LinearNDInterpolator(data_coolant_mu[:, 0:2], data_coolant_mu[:, 2])
    function_k   = LinearNDInterpolator(data_coolant_k[:, 0:2], data_coolant_k[:, 2])
    function_rho = LinearNDInterpolator(data_coolant_rho[:, 0:2], data_coolant_rho[:, 2])

    return function_cp, function_mu, function_k, function_rho


def get_interpolation_function_gas_properties ():

    # function_T        --> [K]         - function to get the temperature of the gas in the nozzle
    # function_cp       --> [kJ/kgK]    - function to get the specific heat capacity of the gas in the nozzle
    # function_k        --> [mW/mK]     - function to get the thermal conductivity of the gas in the nozzle
    # function_mu       --> [uPa*s]     - function to get the viscosity of the gas in the nozzle
    # function_rho      --> [kg/m3]     - function to get the density of the gas in the nozzle
    # function_gamma    --> [-]         - function to get the specific heat ratio of the gas in the nozzle
    # function_c        --> [m/s]       - function to get the speed of sound of the gas in the nozzle
    # function_p        --> [Pa]        - function to get the pressure of the gas in the nozzle

    # Open .csv file using pandas convert to numpy array
    data_gas_T      = pd.read_csv("{}_T.csv".format(CONF_GAS_FILES), header=0).to_numpy()       # expansion_ratio
    data_gas_cp     = pd.read_csv("{}_cp.csv".format(CONF_GAS_FILES), header=0).to_numpy()      # expansion_ratio
    data_gas_mu     = pd.read_csv("{}_mu.csv".format(CONF_GAS_FILES), header=0).to_numpy()      # expansion_ratio
    data_gas_k      = pd.read_csv("{}_k.csv".format(CONF_GAS_FILES), header=0).to_numpy()       # expansion_ratio
    data_gas_rho    = pd.read_csv("{}_rho.csv".format(CONF_GAS_FILES), header=0).to_numpy()     # expansion_ratio
    data_gas_gamma  = pd.read_csv("{}_gamma.csv".format(CONF_GAS_FILES), header=0).to_numpy()   # expansion_ratio
    data_gas_c      = pd.read_csv("{}_c.csv".format(CONF_GAS_FILES), header=0).to_numpy()       # expansion_ratio
    data_gas_p      = pd.read_csv("{}_p.csv".format(CONF_GAS_FILES), header=0).to_numpy()       # expansion_ratio

    # Interolation for the noozle region
    function_T_nozzle      = interp1d(data_gas_T[1:, 0], data_gas_T[1:, 1])
    function_cp_nozzle     = interp1d(data_gas_cp[1:, 0], data_gas_cp[1:, 1])
    function_mu_nozzle     = interp1d(data_gas_mu[1:, 0], data_gas_mu[1:, 1])
    function_k_nozzle      = interp1d(data_gas_k[1:, 0], data_gas_k[1:, 1])
    function_rho_nozzle    = interp1d(data_gas_rho[1:, 0], data_gas_rho[1:, 1])
    function_gamma_nozzle  = interp1d(data_gas_gamma[1:, 0], data_gas_gamma[1:, 1])
    function_c_nozzle      = interp1d(data_gas_c[1:, 0], data_gas_c[1:, 1])
    function_p_nozzle      = interp1d(data_gas_p[1:, 0], data_gas_p[1:, 1])

    functions_nozzle =  [function_T_nozzle, function_cp_nozzle, function_mu_nozzle, function_k_nozzle, function_rho_nozzle, function_gamma_nozzle, function_c_nozzle, function_p_nozzle]

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
    data_gas_p[0,0]         = expansion_ratio_cc

    function_T_cc           = interp1d(data_gas_T[:2, 0], data_gas_T[:2, 1])
    function_cp_cc          = interp1d(data_gas_cp[:2, 0], data_gas_cp[:2, 1])
    function_mu_cc          = interp1d(data_gas_mu[:2, 0], data_gas_mu[:2, 1])
    function_k_cc           = interp1d(data_gas_k[:2, 0], data_gas_k[:2, 1])
    function_rho_cc         = interp1d(data_gas_rho[:2, 0], data_gas_rho[:2, 1])
    function_gamma_cc       = interp1d(data_gas_gamma[:2, 0], data_gas_gamma[:2, 1])
    function_c_cc           = interp1d(data_gas_c[:2, 0], data_gas_c[:2, 1])
    function_p_cc           = interp1d(data_gas_p[:2, 0], data_gas_p[:2, 1])

    functions_cc =  [function_T_cc, function_cp_cc, function_mu_cc, function_k_cc, function_rho_cc, function_gamma_cc, function_c_cc, function_p_cc]

    return functions_cc, functions_nozzle


def graph_interpolation_gas_properties ():

    # Graph the interpolation functions

    # Get the interpolation functions
    functions_cc, functions_nozzle = get_interpolation_function_gas_properties()

    # Create the x-axis
    x_axis = np.linspace(2, 1, 100)

    # Graph the gas properties
    plt.figure()
    plt.plot(x_axis, functions_cc[0](x_axis), label="T_cc")
    plt.plot(x_axis, functions_cc[1](x_axis), label="cp_cc")
    plt.plot(x_axis, functions_cc[2](x_axis), label="mu_cc")
    plt.plot(x_axis, functions_cc[3](x_axis), label="k_cc")
    plt.plot(x_axis, functions_cc[4](x_axis), label="rho_cc")
    plt.plot(x_axis, functions_cc[5](x_axis), label="gamma_cc")
    plt.plot(x_axis, functions_cc[6](x_axis), label="c_cc")
    plt.legend()
    plt.show()

    x_axis = np.linspace(1, 70, 100)

    plt.figure()
    plt.plot(x_axis, functions_nozzle[0](x_axis), label="T_nozzle")
    plt.plot(x_axis, functions_nozzle[1](x_axis), label="cp_nozzle")
    plt.plot(x_axis, functions_nozzle[2](x_axis), label="mu_nozzle")
    plt.plot(x_axis, functions_nozzle[3](x_axis), label="k_nozzle")
    plt.plot(x_axis, functions_nozzle[4](x_axis), label="rho_nozzle")
    plt.plot(x_axis, functions_nozzle[5](x_axis), label="gamma_nozzle")
    plt.plot(x_axis, functions_nozzle[6](x_axis), label="c_nozzle")
    plt.legend()
    plt.show()


if __name__ == "__main__":
    graph_interpolation_gas_properties()