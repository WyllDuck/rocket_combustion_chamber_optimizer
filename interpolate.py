# Tools
from scipy.interpolate import LinearNDInterpolator, interp1d
import pandas as pd

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

    return