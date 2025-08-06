import matplotlib.pyplot as plt
import numpy as np

# Define the logistic function
def logistic(x, k, x0):
    return 1 / (1 + np.exp(-k * (x - x0)))

# Parameters for the logistic function
k = 10  # Steepness of the curve
x0 = 0.5  # Midpoint

# Define the grid
x = np.linspace(0, 1, 100)
y = np.linspace(0, 1, 100)
X, Y = np.meshgrid(x, y)

# Create logistic profiles for X and Y
logistic_X = logistic(X, k, x0)
logistic_Y = logistic(Y, k, x0)

# Sum the logistic profiles for Z
Z = logistic_X + logistic_Y

# Invert Z so that the maximum is at (0, 0) and the minimum is at (1, 1)
Z = 2 - Z  # Subtract Z from 2 to invert it

# Plot
fig = plt.figure(figsize=(8, 6))
ax = fig.add_subplot(111, projection='3d')
ax.plot_surface(X, Y, Z, cmap='viridis')

# Labels
ax.set_xlabel('|ψ1⟩')
ax.set_ylabel('|ψ2⟩')
ax.set_zlabel('E(ψ)')
ax.set_zticks([])

plt.show()
