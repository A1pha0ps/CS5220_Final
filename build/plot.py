import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import subprocess

# Function to parse the data from the file
def parse_data(filename):
    data = []
    with open(filename, 'r') as file:
        matrix = []
        for line in file:
            if line.strip():  # Ignore empty lines
                row = [float(val) for val in line.split()]
                matrix.append(row)
            else:
                data.append(matrix)
                matrix = []
        if matrix:  # Append the last matrix if not empty
            data.append(matrix)
    return data

# subprocess.run("./serial", shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)

# Parse the data from the file
data = parse_data("ex_field.txt")

N = 500 # number of timesteps

# animate self.E_x_t as a .gif file.
# N: number of total steps to save as .gif animation.
E_x_t = np.array(data[-N:])

fig, ax = plt.subplots(figsize = (5, 5))
cax = ax.pcolormesh(np.arange(120), np.arange(120), E_x_t[0].T, 
                    vmin = np.min(E_x_t), vmax = np.max(E_x_t), 
                    shading = "auto", cmap = "bwr")
plt.axis("equal")
plt.colorbar(cax, ax=ax)
plt.grid(True)

def animate(i):
    cax.set_array(E_x_t[i].T.flatten())

anim = FuncAnimation(fig, animate, interval = 50, frames = len(E_x_t) - 1)
anim.save("fdtd_2d_E_animation.gif", writer = "pillow")
# plt.show()
