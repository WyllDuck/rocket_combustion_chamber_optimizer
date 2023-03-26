# Tools
import numpy as np
from stl import mesh
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt

# Global variables
from conf import *


# Using an existing stl file:
your_mesh = mesh.Mesh.from_file('stage_1.stl')
your_mesh.points /= 1000 # convert to meters

# stage 1
your_mesh.points[:,0] = your_mesh.points[:,0] * -1


# print and save in numpy array points in the z-y plane
points = None

for i in range(0, len(your_mesh.points)):

# interpolate spline line
    
    if points is not None:
        points = np.vstack([points, np.array([your_mesh.points[i][0], your_mesh.points[i][1]])])
    else:
        points = np.array([your_mesh.points[i][0], your_mesh.points[i][1]])


final_points = None

x = np.unique(points[:,0])
for x_ in x:

    location    = np.where(points[:,0] == x_)[0]
    y_max       = max(points[location][:,1])

    if x_ == 0:
        print ("throat D: {}".format(y_max * 2))


    if y_max < 0:
        continue

    if final_points is not None:
        final_points = np.vstack([final_points, np.array([x_, y_max])])
    else:
        final_points = np.array(np.array([x_, y_max]))

f = interp1d(final_points[:,0], final_points[:,1], kind='linear')

"""
COMPILE DATA
"""


x_min = min(final_points[:,0])
x_max = max(final_points[:,0])

x_int1 = np.linspace(x_min, 0, 100)
x_int2 = np.linspace(0, x_max, 100)[1:] # remove the first point because it is the same as the last point of the first interval

data_final = np.zeros([len(x_int1) + len(x_int2), 3])

data_final[len(x_int1):, 0] = x_int2
data_final[len(x_int1):, 1] = f(x_int2)
data_final[len(x_int1):, 2] = 2

data_final[:len(x_int1), 0] = x_int1
data_final[:len(x_int1), 1] = f(x_int1)
data_final[:len(x_int1), 2] = 1

# Save the data into a file 
data_final[:, 1] = data_final[:, 1] * 2     # Radius is doubled because we need the diameter
np.savetxt(GEOMETRY_FILE, data_final[::-1], delimiter=',')

loc_points_cc = np.where(data_final[:,2] == 1)[0]
loc_points_no = np.where(data_final[:,2] == 2)[0]

# plot the points
plt.plot(data_final[loc_points_cc,0], data_final[loc_points_cc,1], 'o')
plt.plot(data_final[loc_points_no,0], data_final[loc_points_no,1], '-')

plt.plot(final_points[:,0], final_points[:,1] * 2, 'o')

plt.show()