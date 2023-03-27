# Other files
from solver import Solver
from plot import *

# Global variables
from conf import *


""" MAIN """
def main ():

    dir_ = "img"

    solver = Solver()
    solver.solve()

    plot_thickness(solver, save=True, filename="{}/thickness.svg".format(dir_))
    plot_temperature(solver, save=False, filename="{}/temperature.svg".format(dir_))
    plot_pressure(solver, save=False, filename="{}/pressure.svg".format(dir_))
    plot_velocity(solver, save=False, filename="{}/velocity.svg".format(dir_))
    plot_heat_flux(solver, save=False, filename="{}/heat_flux.svg".format(dir_))
    plot_heat_flux_cc(solver, save=False, filename="{}/heat_conv_coeffients_cc.svg".format(dir_))
    

if __name__ == "__main__":
    main()
