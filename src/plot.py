# Tools
import numpy as np
import matplotlib.pyplot as plt

# global variables
from conf import *


# Auxiliary functions
def get_xy_coord ():

    # Open the geometry file
    data = np.loadtxt(GEOMETRY_FILE, delimiter=',')
    xy = np.zeros([len(data) - 1, 2])

    for i in range(len(data) - 1):
        xy[i, 0] = 0.5*(data[i + 1, 0] + data[i, 0])
        xy[i, 1] = 0.5*(data[i + 1, 1] + data[i, 1])

    return xy


# plot temperature using matplotlib.pyplot and numpy of the coolant circuit given a Solver object
def plot_temperature (solver, save=False, filename="img/temperature.svg"):

    xy = get_xy_coord()
    x_ = xy[:,0]
    y_ = xy[:,1]

    T = np.zeros([len(x_), 3])

    for i in range(len(x_)):
        T[i, 0] = solver.sections[i].T_in
        T[i, 1] = solver.sections[i].T_wi_
        T[i, 2] = solver.sections[i].T_wo_

    plt.plot(x_, T[:, 0], label="T coolant")
    plt.plot(x_, T[:, 1], label="T coolant wall")
    plt.plot(x_, T[:, 2], label="T hot gases wall")

    # plot horizontal line at MELTING_T
    plt.axhline(y=MELTING_T, color='r', linestyle='--', label="Melting Temperature")

    # BACKGROUND NOZZLE for reference
    ylim_min, ylim_max = plt.gca().get_ylim()
    y_ = (y_ - min(y_)) / (max(y_) - min(y_))   # normalize y_ to [0, 1]
    y_ = y_ * (ylim_max - ylim_min) + ylim_min  # scale y_ to [ylim_min, ylim_max]

    plt.plot(x_, y_, color="gray", linestyle="--", label="Nozzle Profile")

    plt.gca().set_xlabel("x [m]")
    plt.gca().set_ylabel("T [K]")

    plt.legend()
    plt.grid()

    if save:
        plt.savefig(filename)
    plt.show()


def plot_pressure (solver, save=False, filename="img/pressure.svg"):

    xy = get_xy_coord()
    x_ = xy[:,0]
    y_ = xy[:,1]

    p = np.zeros(len(x_))

    for i in range(len(x_)):
        p[i] = solver.sections[i].p_out

    plt.plot(x_, p, label="Pressure Coolant")
    
    # BACKGROUND NOZZLE for reference
    ylim_min, ylim_max = plt.gca().get_ylim()
    y_ = (y_ - min(y_)) / (max(y_) - min(y_))   # normalize y_ to [0, 1]
    y_ = y_ * (ylim_max - ylim_min) + ylim_min  # scale y_ to [ylim_min, ylim_max]

    plt.plot(x_, y_, color="gray", linestyle="--", label="Nozzle Profile")

    plt.gca().set_xlabel("x [m]")
    plt.gca().set_ylabel("P [bar]")

    plt.legend()
    plt.grid()

    if save:
        plt.savefig(filename)
    plt.show()


def plot_velocity (solver, save=False, filename="img/velocity.svg"):

    xy = get_xy_coord()
    x_ = xy[:,0]
    y_ = xy[:,1]

    v = np.zeros(len(x_))

    for i in range(len(x_)):
        v[i] = solver.sections[i].v

    plt.plot(x_, v, label="Coolant Velocity")
    
    # BACKGROUND NOZZLE for reference
    ylim_min, ylim_max = plt.gca().get_ylim()
    y_ = (y_ - min(y_)) / (max(y_) - min(y_))   # normalize y_ to [0, 1]
    y_ = y_ * (ylim_max - ylim_min) + ylim_min  # scale y_ to [ylim_min, ylim_max]

    plt.plot(x_, y_, color="gray", linestyle="--", label="Nozzle Profile")

    plt.gca().set_xlabel("x [m]")
    plt.gca().set_ylabel("v [m/s]")

    plt.legend()
    plt.grid()

    if save:
        plt.savefig(filename)
    plt.show()


def plot_heat_flux (solver, save=False, filename="img/heat_flux.svg"):
    
    xy = get_xy_coord()
    x_ = xy[:,0]
    y_ = xy[:,1]

    Q = np.zeros(len(x_))

    for i in range(len(x_)):
        Q[i] = solver.sections[i].Q

    plt.plot(x_, Q, label="Heat Flus Q")
    
    # BACKGROUND NOZZLE for reference
    ylim_min, ylim_max = plt.gca().get_ylim()
    y_ = (y_ - min(y_)) / (max(y_) - min(y_))   # normalize y_ to [0, 1]
    y_ = y_ * (ylim_max - ylim_min) + ylim_min  # scale y_ to [ylim_min, ylim_max]

    plt.plot(x_, y_, color="gray", linestyle="--", label="Nozzle Profile")

    plt.gca().set_xlabel("x [m]")
    plt.gca().set_ylabel("Q [W/m^2]")

    plt.legend()
    plt.grid()

    if save:
        plt.savefig(filename)
    plt.show()


def plot_heat_flux_cc (solver, save=False, filename="img/heat_conv_coeffients_cc.svg"):
    
    xy = get_xy_coord()
    x_ = xy[:,0]
    y_ = xy[:,1]

    h = np.zeros([len(x_), 2])

    for i in range(len(x_)):
        h[i, 0] = solver.sections[i].h_cc
        h[i, 1] = solver.sections[i].h_co

    plt.plot(x_, h[:, 0], label="heat convection coeff. combustion chamber")
    plt.plot(x_, h[:, 1], label="heat convection coeff. cooling channel")
    
    # BACKGROUND NOZZLE for reference
    ylim_min, ylim_max = plt.gca().get_ylim()
    y_ = (y_ - min(y_)) / (max(y_) - min(y_))   # normalize y_ to [0, 1]
    y_ = y_ * (ylim_max - ylim_min) + ylim_min  # scale y_ to [ylim_min, ylim_max]

    plt.plot(x_, y_, color="gray", linestyle="--", label="nozzle profile")

    plt.gca().set_xlabel("x [m]")
    plt.gca().set_ylabel("h [W/m^2/K]")

    plt.legend()
    plt.grid()

    if save:
        plt.savefig(filename)
    plt.show()
