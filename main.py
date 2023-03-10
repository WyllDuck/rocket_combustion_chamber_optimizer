# Other files
from solver import Solver
from plot import *

# Global variables
from conf import *


""" MAIN """
def main ():

    solver = Solver()
    solver.solve()

    plot_temperature(solver)
    plot_pressure(solver)
    plot_velocity(solver)
    plot_heat_flux(solver)
    plot_heat_flux_cc(solver)
    

if __name__ == "__main__":
    main()
