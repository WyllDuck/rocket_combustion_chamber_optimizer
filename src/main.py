# Other files
from solver import Solver
from plot import *

# Global variables
from conf import *


""" MAIN """
def main ():

    solver = Solver()
    solver.solve()

    plot_temperature(solver, save=True)
    plot_pressure(solver, save=True)
    plot_velocity(solver, save=True)
    plot_heat_flux(solver, save=True)
    plot_heat_flux_cc(solver, save=True)
    

if __name__ == "__main__":
    main()
