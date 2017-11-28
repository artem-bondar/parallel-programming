from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import numpy as np
 
fig = plt.figure(figsize=(10,6))
ax = fig.add_subplot(111, projection='3d')
ax.set_xlabel("Threads")
ax.set_ylabel("Size")
ax.set_zlabel("Time")

x_data = []
y_data = []
z_data = []
 
file = open("Regatta OpenMP.txt","r")
results = file.readlines()

for line in results:
    data = line.split(" ")
    x_data.append(int(data[0]))
    y_data.append(int(data[1]))
    z_data.append(int(data[2].split(".")[0]))

dx = [0.2 for i in z_data]
dy = [100 for i in z_data]

ax.bar3d( x_data, y_data, np.zeros(len(z_data)), dx, dy, z_data, color='#0099cc')

plt.show()