import numpy as np
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt

# Global variables
from conf import *


# Function that opens file and return x, y coordinates in a numpy array
def get_data(filename):
    data = np.loadtxt(filename, delimiter='\t')
    return data[:, 1:3]


def main ():

    """
    COMBUSTION CHAMBER
    """
    x_int1 = np.linspace(0, LENGHT_CC, 100)

    """
    NOZZLE EXTENSION
    """
    data = get_data("geo/2ndStage_eps150_Ma4_678_gamma_1_1982_contour.txt") * DI_TH / 2 # Dimensionalize the data

    CUT_NOZZLE_R = 0.8 #0.25
    for i in range(len(data)):
        if (data[i, 1] - CUT_NOZZLE_R) > 0:
            break
    data = data[:i+1, :]

    data[:, 0] = data[:, 0] + LENGHT_CC # move x axis to the end of the combustion chamber

    # Interpolate the data
    f = interp1d(data[:, 0], data[:, 1], kind='cubic')
    x_int2 = np.linspace(data[0, 0], data[-1, 0], 100)

    """
    COMPILE DATA
    """
    data_final = np.zeros([len(x_int1) + len(x_int2), 3])
    
    data_final[len(x_int1):, 0] = x_int2
    data_final[len(x_int1):, 1] = f(x_int2)
    data_final[len(x_int1):, 2] = 2

    data_final[:len(x_int1), 0] = x_int1
    data_final[:len(x_int1), 1] = DI_CC / 2 # Radius of the combustion chamber
    data_final[:len(x_int1), 2] = 1

    """
    PLOT DATA
    """
    # plot using othorgonal axes
    plt.axis('equal')
    plt.plot(data_final[:, 0], data_final[:, 1], color='black')
    plt.plot(data_final[:, 0], -data_final[:, 1], color='black')
    plt.show()

    # Save the data into a file 
    data_final[:, 1] = data_final[:, 1] * 2     # Radius is doubled because we need the diameter
    np.savetxt(GEOMETRY_FILE, data_final[::-1], delimiter=',')


if __name__ == "__main__":
    main()