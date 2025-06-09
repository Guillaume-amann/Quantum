import pulp
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import random

# Bin dimensions
W, H = 4, 3

# Items: (width, height, count)
items = [
    (1, 1, 3),  # 3 small squares
    (2, 1, 4),  # 5 horizontal rectangles
]

reward = 10
costs = [5, 3]  # cost for 1x1 is 5, for 2x1 is 3
penalty = -10

model = pulp.LpProblem("2D_Rectangular_Bin_Packing", pulp.LpMaximize)

# Decision vars: x[item_id][i][j] = 1 if item placed at (i,j)
x = {}
for k, (iw, ih, count) in enumerate(items):
    for i in range(H - ih + 1):
        for j in range(W - iw + 1):
            x[k, i, j] = pulp.LpVariable(f"x_{k}_{i}_{j}", cat="Binary")

# Cell coverage matrix: cell_used[i][j] = 1 if any item covers (i,j)
cell_used = {}
for i in range(H):
    for j in range(W):
        cell_used[i, j] = pulp.LpVariable(f"used_{i}_{j}", cat="Binary")

# Constraint: each cell used at most once
for i in range(H):
    for j in range(W):
        covering_items = []
        for k, (iw, ih, _) in enumerate(items):
            for di in range(max(0, i - ih + 1), min(i + 1, H - ih + 1)):
                for dj in range(max(0, j - iw + 1), min(j + 1, W - iw + 1)):
                    if (i - di < ih) and (j - dj < iw):
                        covering_items.append(x[k, di, dj])
        model += cell_used[i, j] == pulp.lpSum(covering_items)

# Constraint: item usage limited to their count
for k, (_, _, count) in enumerate(items):
    model += pulp.lpSum([x[k, i, j] for (k2, i, j) in x if k2 == k]) <= count

# Objective: maximize reward - cost - penalties
filled_value = reward
empty_value = penalty
objective = pulp.lpSum([
    cell_used[i, j] * reward + (1 - cell_used[i, j]) * penalty
    for i in range(H) for j in range(W)
]) - pulp.lpSum([
    x[k, i, j] * costs[k]
    for (k, i, j) in x
])

model += objective

# Solve
model.solve()
print("Status:", pulp.LpStatus[model.status])
print("Objective value:", pulp.value(model.objective))



# Build grid
grid = np.full((H, W), '.', dtype=str)

for (k, i, j), var in x.items():
    if pulp.value(var) == 1:
        w, h, _ = items[k]
        for di in range(h):
            for dj in range(w):
                grid[i + di][j + dj] = chr(65 + k)  # A, B, etc.

print("Grid:")
print(grid)

# Plot
def plot_grid(grid):
    H, W = grid.shape
    fig, ax = plt.subplots()

    # Identify unique item instances (exclude empty cells)
    unique_items = sorted(set(grid.flatten()) - set('.'))

    # Map each unique item to a unique index and random color
    label_to_int = {label: idx+1 for idx, label in enumerate(unique_items)}
    color_grid = np.zeros((H, W))

    for i in range(H):
        for j in range(W):
            if grid[i][j] != '.':
                color_grid[i][j] = label_to_int[grid[i][j]]

    # Generate a random color for each unique item instance
    cmap_colors = ["white"]  # index 0 for empty cell
    for _ in unique_items:
        color = [random.random() for _ in range(3)]
        cmap_colors.append(mcolors.to_hex(color))

    cmap = mcolors.ListedColormap(cmap_colors)

    ax.imshow(color_grid, cmap=cmap, origin="upper")

    ax.set_xticks(np.arange(-0.5, W, 1), minor=True)
    ax.set_yticks(np.arange(-0.5, H, 1), minor=True)
    ax.grid(which="minor", color="black", linewidth=1)
    ax.set_xticks([])
    ax.set_yticks([])

    # Draw item labels on each cell
    for i in range(H):
        for j in range(W):
            if grid[i][j] != '.':
                ax.text(j, i, grid[i][j], ha="center", va="center", color="black", fontsize=12, fontweight='bold')

    plt.title("2D Bin Packing — Random Color per Item Instance")
    plt.show()

plot_grid(grid)