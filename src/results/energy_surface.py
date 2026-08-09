import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

df = pd.read_csv("energy_surface.csv")
fig = plt.figure(figsize=(10, 6))
ax = fig.add_subplot(111, projection='3d')
ax.plot_trisurf(df['gamma'], df['alpha'], df['energy'], cmap='viridis')
ax.set_xlabel('Gamma')
ax.set_ylabel('Alpha')
ax.set_zlabel('Energy')
ax.set_title('Energy Surface E(γ, α)')
plt.tight_layout()
plt.show()