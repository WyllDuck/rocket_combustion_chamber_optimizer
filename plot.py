# Tools
import numpy as np
import matplotlib.pyplot as plt

# global variables
from conf import *


# Auxiliary functions
def get_x_coord ():

    # Open the geometry file
    data = np.loadtxt(GEOMETRY_FILE, delimiter=',')
    x = np.zeros(len(data) - 1)

    for i in range(len(data) - 1):
        x[i] = 0.5*(data[i + 1, 0] + data[i, 0])

    return x


# plot temperature using matplotlib.pyplot and numpy of the coolant circuit given a Solver object
def plot_temperature (solver, save=False, filename="img/temperature.png"):

    x_ = get_x_coord()
    T = np.zeros([len(x_), 3])

    for i in range(len(x_)):
        T[i, 0] = solver.sections[i].T_in
        T[i, 1] = solver.sections[i].T_wi_
        T[i, 2] = solver.sections[i].T_wo_

    plt.plot(x_, T[:, 0], label="T_in")
    plt.plot(x_, T[:, 1], label="T_wi_")
    plt.plot(x_, T[:, 2], label="T_wo_")

    plt.legend()
    if save:
        plt.savefig(filename)
    plt.show()

def plot_pressure (solver, save=False, filename="img/pressure.png"):

    x_ = get_x_coord()
    p = np.zeros(len(x_))

    for i in range(len(x_)):
        p[i] = solver.sections[i].p_out

    plt.plot(x_, p, label="p_out")
    
    plt.legend()
    if save:
        plt.savefig(filename)
    plt.show()

def plot_velocity (solver, save=False, filename="img/velocity.png"):

    x_ = get_x_coord()
    v = np.zeros(len(x_))

    for i in range(len(x_)):
        v[i] = solver.sections[i].v

    plt.plot(x_, v, label="v")
    
    plt.legend()
    if save:
        plt.savefig(filename)
    plt.show()


def plot_heat_flux (solver, save=False, filename="img/heat_flux.png"):
    
    x_ = get_x_coord()
    Q = np.zeros(len(x_))

    for i in range(len(x_)):
        Q[i] = solver.sections[i].Q

    plt.plot(x_, Q, label="Q")
    
    plt.legend()
    if save:    
        plt.savefig(filename)
    plt.show()


def plot_heat_flux_cc (solver, save=False, filename="img/heat_flux_cc.png"):
    
    x_ = get_x_coord()
    h = np.zeros([len(x_), 2])

    for i in range(len(x_)):
        h[i, 0] = solver.sections[i].h_cc
        h[i, 1] = solver.sections[i].h_co

    plt.plot(x_, h[:, 0], label="h_cc")
    plt.plot(x_, h[:, 1], label="h_co")
    
    plt.legend()
    if save:
        plt.savefig(filename)
    plt.show()
