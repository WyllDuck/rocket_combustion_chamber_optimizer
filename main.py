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
    

if __name__ == "__main__":
    main()
