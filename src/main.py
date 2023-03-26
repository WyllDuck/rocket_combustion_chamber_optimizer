# Other files
from solver import Solver
from plot import *

# Global variables
from conf import *


""" MAIN """
def main ():

    solver = Solver()
    solver.solve()

    plot_temperature(solver, save=False)
    plot_pressure(solver, save=False)
    plot_velocity(solver, save=False)
    plot_heat_flux(solver, save=False)
    plot_heat_flux_cc(solver, save=False)
    

if __name__ == "__main__":
    main()
