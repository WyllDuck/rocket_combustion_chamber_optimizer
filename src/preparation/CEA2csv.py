# TOOLS
import numpy as np

# Global variables
from conf import *


INDEXES = ( [99, "Ae/At", "expansion_ratio", "-"],

            [64, "P, BAR", "p", "bar"],
            [65, "T, K", "T", "K"],
            [66, "RHO, KG/CU M", "rho", "kg/m^3"],
            [75, "Cp, KJ/(KG)(K)", "cp", "kJ/(kg K)"],
            [76, "GAMMAs", "gamma", "-"],
            [77, "SON VEL,M/SEC", "c", "m/s"],
            [83, "VISC,MILLIPOISE", "mu", "millipoise"],


            # WITH EQUILIBRIUM REACTIONS
            #[87, "Cp, KJ/(KG)(K)"],
            [88, "CONDUCTIVITY", "k", "mW/(cm K)"],


            # WITH FROZEN REACTIONS
            #[93, "Cp, KJ/(KG)(K)"],
            #[94, "CONDUCTIVITY", "k", "mW/(cm K)"],
            )





# Read data from CEA file
def get_data_from_CEA (cea_file, all_content = None):

    # if all_content is None, create a new list
    # else, use the list that was passed and expand data
    if all_content is None:
        all_content = [None]*(len(INDEXES) - 1)

    # Open the file
    f = open(cea_file, "r")
    lines = f.readlines()

    # Get header
    header  = process_line(lines, 0)
    header  = np.concatenate([np.array([0]), header])   # add 0 to the beginning of the array to signal values inside CC

    for i in range (len(all_content)):
        cont    = process_line(lines, i+1)

        # create content matrix
        data    = np.zeros((len(header), 2))
        data[:,0]   = header
        data[:,1]   = cont

        if not np.any(all_content[i]):
            all_content[i] = data
        else:
            all_content[i] = np.vstack([all_content[i], data])

    return all_content


# Save data to file CSV files
def save_data (file_name, content):

    # save data to file with name
    for i in range(1, len(INDEXES)):
        file_name_  = file_name + "_" + INDEXES[i][2] + ".csv"
        header_     = "{}[{}]".format(INDEXES[0][2], INDEXES[0][3]) + ", " + "{}[{}]".format(INDEXES[i][2], INDEXES[i][3])

        np.savetxt(file_name_, content[i-1], delimiter=",", header=header_)

        
# Process a line from CEA file
def process_line (data, i):

    idx, header, _, _     = INDEXES[i]
    print("idx: {}, header: {}".format(idx, header))
    idx     -= 1 # because python starts at 0

    line    = data[idx]
    line    = line.replace("\n", "")
    line    = line.replace(header, "")

    line    = line.split()

    line_ = []
    for i in range(len(line)):
        if not "+" in line[i] and not "-" in line[i]:
            line_.append(float(line[i]))
        else:
            idx_sign = line[i].find("+")
            if idx_sign == -1:
                idx_sign = line[i].find("-")
            line_.append(float(line[i][:idx_sign] + "e" + line[i][idx_sign:]))

    print(line_)
    return np.array(line_)


# Main
if __name__ == "__main__":
    
    content = None
    for i in range (len(CEA_FILES)):
        content = get_data_from_CEA(CEA_FILES[i], content)

    # save file 
    save_data(CONF_GAS_FILES, content)