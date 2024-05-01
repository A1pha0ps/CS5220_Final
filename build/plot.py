import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from matplotlib.patches import Rectangle

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
data = parse_data("e_field_before.txt")

N = 200 # number of timesteps

# animate self.E_x_t as a .gif file.
# N: number of total steps to save as .gif animation.
E_x_t = np.array(data[-N:])

fig, ax = plt.subplots(figsize = (5, 5))
cax = ax.pcolormesh(np.arange(48), np.arange(48), E_x_t[0].T, 
                    vmin = np.min(E_x_t), vmax = np.max(E_x_t), 
                    shading = "auto", cmap = "bwr")
plt.axis("equal")
plt.colorbar(cax, ax=ax)
plt.grid(True)

# Define rectangle parameters
# rect_x = 75
# rect_y = 0
# rect_width = 150
# rect_height = 300
# rect_edgecolor = 'black'
# rect_facecolor = 'none'

# Add rectangle patch
# rect = Rectangle((rect_x, rect_y), rect_width, rect_height, linewidth=2, edgecolor=rect_edgecolor, facecolor=rect_facecolor)
# ax.add_patch(rect)

def animate(i):
    cax.set_array(E_x_t[i].T.flatten())

anim = FuncAnimation(fig, animate, interval = 50, frames = len(E_x_t) - 1)
anim.save("E_before.gif", writer = "pillow")
# plt.show()
