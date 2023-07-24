# -*- coding: utf-8 -*-
"""Created on Wed Jun 14 14:36:20 2023.

@author: dpm42
"""

import seaborn as sns

# %%
sns.set(font_scale=1.25, style='whitegrid')
import matplotlib.pyplot as plt
import numpy as np


def calculate_tortuosity(x, y):
    total_distance = np.sum(np.sqrt(np.diff(x) ** 2 + np.diff(y) ** 2))
    straight_distance = np.sqrt((x[-1] - x[0]) ** 2 + (y[-1] - y[0]) ** 2)
    tortuosity = total_distance / straight_distance
    return tortuosity


def plot_tortuosity(non_directness, ax):
    num_points = 100
    x = np.zeros(num_points)
    y = np.linspace(0, 10, num_points)

    # Adding non-directness to the line
    deviation = 0.2 * non_directness
    x += np.random.uniform(-deviation, deviation, size=num_points)

    tortuosity = calculate_tortuosity(x, y)

    ax.plot(y, x, 'k', linewidth=2)
    ax.set_ylabel('')
    ax.set_title(f"Tortuosity: {tortuosity:.2f}", pad=-17, y=1.001)
    ax.set_ylim(-1, 1)  # Set the same x limit for all plots
    ax.set_aspect('equal')
    ax.grid(False)


# Generating plots for different levels of non-directness
non_directness_values = [0, 0.1, 0.2]
# Create subplots
fig, axs = plt.subplots(len(non_directness_values), 1, figsize=(6, 6), sharex=True)
# axs[0].set_ylabel('Z (a.u.)')
for i, non_directness in enumerate(non_directness_values):
    plot_tortuosity(non_directness, axs[i])
plt.subplots_adjust(hspace=-0.8)
plt.tight_layout()
plt.show()
