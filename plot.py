# Tools
import numpy as np
import matplotlib.pyplot as plt

# global variables
from conf import *


# plot temperature using matplotlib.pyplot and numpy of the coolant circuit given a Solver object
def plot_temperature (solver):

    T = np.zeros([N_SECTIONS, 3])
    dx_ = np.linspace(0, LENGHT_CC, N_SECTIONS)

    for i in range(N_SECTIONS):
        T[i, 0] = solver.sections[i].T_in
        T[i, 1] = solver.sections[i].T_wi_
        T[i, 2] = solver.sections[i].T_wo_

    plt.plot(dx_, T[:, 0], label="T_in")
    plt.plot(dx_, T[:, 1], label="T_wi_")
    plt.plot(dx_, T[:, 2], label="T_wo_")

    plt.legend()
    plt.show()

def plot_pressure (solver):

    p = np.zeros(N_SECTIONS)
    dx_ = np.linspace(0, LENGHT_CC, N_SECTIONS)

    for i in range(N_SECTIONS):
        p[i] = solver.sections[i].p_out

    plt.plot(dx_, p, label="p_out")
    
    plt.legend()
    plt.show()

def plot_velocity (solver):

    v = np.zeros(N_SECTIONS)
    dx_ = np.linspace(0, LENGHT_CC, N_SECTIONS)

    for i in range(N_SECTIONS):
        v[i] = solver.sections[i].v

    plt.plot(dx_, v, label="v")
    
    plt.legend()
    plt.show()
