# Tools
import pyromat as pm
import yaml


# FIXED GLOBAL PARAMETERS
INLET_T     = 200   # [K]
INLET_p     = 130   # [bar]

THERMAL_K   = 385   # [W/mK]
CC_LENGTH   = 0.3   # [m]
N_SECTIONS  = 20    # [º]

N_CHANNELS  = 20    # [º]


class CombustionChamber (object):

    def __init__ (self) -> None:
        
        Section.dx = CC_LENGTH / N_SECTIONS

        # Create All cross-sections of the combustion chamber
        self.section = list(range(N_SECTIONS))
        self.section[0] = Section() 

        for i in range (1, N_SECTIONS):
            self.section[i] = Section(self.section[i-1])


    def execute (self):

        for section in self.sections:
            section.get_exit_conditions()

        return 



# Cross section of the combustion chamber 
class Section (object):

    dx = None # value set in CombutionChamber initialization 

    def __init__ (self, previous_section = None) -> None:
        
        self.pre_section = previous_section     # section before this one.
        
        # NOTE: If a previous section is provided then override the values, if not use CC_INLET conditions
        if not self.pre_section:        
            self.T_in   = self.pre_section.T_out
            self.p_in   = self.pre_section.p_out

        # INLET conditions
        else:
            self.T_in   = INLET_T # global variable
            self.p_in   = INLET_p # global variable

        self.T_out = 0
        self.p_out = 0

        # Dimension parameters




    def get_exit_condition (self):

        

    
    def get_pressure_drop ()
        



# Main
def main ():
    combustion_chamber = CombustionChamber()
    combustion_chamber.execute()


if __name__ == "__main__":
    main()
