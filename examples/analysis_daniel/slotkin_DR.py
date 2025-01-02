import matplotlib.pyplot as plt
import numpy as np

# Define the dose levels and response values for 2DEM and 3DM
dose_levels = np.linspace(0, 10, 100)
response_2DEM = 1 / (1 + np.exp(-dose_levels + 4))
response_3DM = 1 / (1 + np.exp(-dose_levels + 4.5))

# Calculate and plot the difference between the two curves (original order)
response_difference = response_3DM - response_2DEM

# Swap the curves
response_2DEM, response_3DM = response_3DM, response_2DEM

# Calculate and plot the difference between the swapped curves
response_difference_swapped = response_3DM - response_2DEM

# Create a square figure
fig, axes = plt.subplots(2, 2, figsize=(8, 8), sharey='col')

# Plot the original curves and response difference
axes[0, 0].plot(dose_levels, response_2DEM, color='orange', label='2DEM')
axes[0, 0].plot(dose_levels, response_3DM, color='purple', label='3DM')
axes[0, 0].set_ylabel('Prop. fibers activated')
axes[0, 0].set_title('Dose-Response Curves')
quants = np.arange(0.05, 1, 0.1)
for q in quants:
    axes[0, 0].axhline(y=q, color='red', linestyle='--', label=f'{"quantile" if q==0.05 else "_"}', alpha=0.4)

axes[0, 0].legend()

axes[0, 1].plot(dose_levels, response_difference, color='blue', label='Difference (3DM - 2DEM)')
axes[0, 1].set_ylabel('Response Difference')
axes[0, 1].set_title('Response Difference')
plt.xticks([0, 10], ['0.05', '0.95'])
# Define the dose levels and response values for 2DEM and 3DM
dose_levels = np.linspace(0, 10, 100)
response_2DEM = 1 / (1 + np.exp(-dose_levels + 3))
response_3DM = 1 / (1 + np.exp(-dose_levels + 7))

# Calculate and plot the difference between the two curves (original order)
response_difference = response_3DM - response_2DEM

# Swap the curves
response_2DEM, response_3DM = response_3DM, response_2DEM

quants = np.arange(0.05, 1, 0.1)
for q in quants:
    axes[1, 0].axhline(y=q, color='red', linestyle='--', label=f'{"quantile" if q==0.05 else "_"}', alpha=0.4)
# Plot the swapped curves and response difference
axes[1, 0].plot(dose_levels, response_2DEM, color='orange', label='2DEM')
axes[1, 0].plot(dose_levels, response_3DM, color='purple', label='3DM')
axes[1, 0].set_xlabel('Threshold')
axes[1, 0].set_ylabel('Prop. fibers activated')

axes[1, 1].plot(dose_levels, response_difference, color='blue', label='Difference (3DM - 2DEM, Swapped)')
axes[1, 1].set_xlabel('Quantile')
axes[1, 1].set_ylabel('Response Difference')
plt.xticks([0, 10], ['0.05', '0.95'])
axes[0, 1].axhline(y=0, color='black', linestyle='--')
axes[1, 1].axhline(y=0, color='black', linestyle='--')

# Adjust layout and show the figure
plt.tight_layout()
plt.show()
# Create a square figure
fig, ax = plt.subplots(1, 1, figsize=(4, 4), sharey='col')
dose_levels = np.linspace(0, 10, 100)

# Plot the original curves and response difference
response_2DEM = 1 / (1 + np.exp(-dose_levels + 5))
ax.plot(dose_levels, response_2DEM, color='black', label='3DM')
ax.set_ylabel('Prop. fibers activated')
ax.set_xlabel('Threshold (mA)')
plt.xticks([0, 10], ["1.2", "3.2"])
