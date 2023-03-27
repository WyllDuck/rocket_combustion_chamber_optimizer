# Tools
import numpy as np
from math import pow, sqrt
from stl import mesh
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt

import sys
import os

dir_path = os.path.abspath(os.path.join(os.path.dirname( __file__ ), '..'))
sys.path.append(dir_path)

# Global variables
from conf import *


# Using an existing stl file:
your_mesh = mesh.Mesh.from_file('stl\stage_2.stl')
your_mesh.points /= 1000 # convert to meters

# stage 1
#your_mesh.points[:,0] = your_mesh.points[:,0] * -1

# print and save in numpy array points in the z-y plane
points = np.array([your_mesh.points[0][0], your_mesh.points[0][1]])
for i in range(1, len(your_mesh.points)):
    d = sqrt(pow(your_mesh.points[i][1], 2) + pow(your_mesh.points[i][2], 2))
    points = np.vstack([points, np.array([your_mesh.points[i][0], d])])

plt.plot(points[:,0], points[:,1], 'o')
plt.axis('equal')
plt.show()

final_points = None

x = np.unique(points[:,0])
for x_ in x:

    location    = np.where(points[:,0] == x_)[0]
    y_max       = max(points[location][:,1])

    if y_max < 0:
        continue

    if final_points is not None:
        final_points = np.vstack([final_points, np.array([x_, y_max])])
    else:
        final_points = np.array(np.array([x_, y_max]))

# Solve error in the combustion chamber profile
x_min   = min(final_points[:,0])
x_min2  = min(final_points[1:,0]) # points to curvature initial point at CC
x_max   = max(final_points[:,0])

extra_points = np.zeros([40, 2])
extra_points[:,0] = np.linspace(x_min, x_min2, 40)
extra_points[:,1] = final_points[0,1]
final_points = np.concatenate([extra_points, final_points[2:]], axis=0) # remove the next 2 points because they look are repeated points in the extra_points array 

f = interp1d(final_points[:,0], final_points[:,1], kind='cubic')

"""
COMPILE DATA
"""




x_int1 = np.linspace(x_min, 0, 1000)
x_int2 = np.linspace(0, x_max, 4000)[1:] # remove the first point because it is the same as the last point of the first interval

data_final = np.zeros([len(x_int1) + len(x_int2), 3])

data_final[len(x_int1):, 0] = x_int2
data_final[len(x_int1):, 1] = f(x_int2)
data_final[len(x_int1):, 2] = 2

data_final[:len(x_int1), 0] = x_int1
data_final[:len(x_int1), 1] = f(x_int1)
data_final[:len(x_int1), 2] = 1

CUT_NOZZLE_R = DI_TH + 0.2 # 6.2 bar line
for i in range(len(x_int1), len(data_final)):
    if (data_final[i, 1] - CUT_NOZZLE_R) > 0 and data_final[i, 2] == 2:
        break
data_final = data_final[:i+1, :]

# Save the data into a file 
data_final[:, 1] = data_final[:, 1] * 2     # Radius is doubled because we need the diameter
np.savetxt(GEOMETRY_FILE, data_final[::-1], delimiter=',')

loc_points_cc = np.where(data_final[:,2] == 1)[0]
loc_points_no = np.where(data_final[:,2] == 2)[0]

# plot the points
plt.plot(data_final[loc_points_cc,0], data_final[loc_points_cc,1], 'r-')
plt.plot(data_final[loc_points_no,0], data_final[loc_points_no,1], 'b-')

plt.plot(final_points[:,0], final_points[:,1] * 2, 'go')

print("DI_CC = {}".format(2*f(x_min)))
print("DI_TH = {}".format(2*f(0)))

plt.show()