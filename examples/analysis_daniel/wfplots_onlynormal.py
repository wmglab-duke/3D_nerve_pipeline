# %%
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import tol_colors as tc

mpl.rcParams['figure.dpi'] = 400
sns.set(style='ticks', context='paper')

# Define the waveform time and amplitude details
t1 = np.linspace(0, 200e-6, 100)  # Active phase for asymmetric waveform
t2 = np.linspace(200e-6, 2.2e-3, 100)  # Rest phase for asymmetric waveform
t2_symmetric = np.linspace(200e-6, 400e-6, 100)  # Time for symmetric waveform

# Asymmetric biphasic square wave
wave1_phase1 = np.ones_like(t1)
wave1_phase2 = -1 / 10 * np.ones_like(t2)  # Scaled to match the area under the curve
wave1 = np.concatenate([wave1_phase1, wave1_phase2])
t_wave1 = np.concatenate([t1, t2])

# Symmetric biphasic waveform
wave2_phase1 = np.ones_like(t1)
wave2_phase2 = -np.ones_like(t1)
wave2 = np.concatenate([wave2_phase1, wave2_phase2])
t_wave2 = np.concatenate([t1, t2_symmetric])

# Adjust waveforms to connect to y=0 at the beginning and end
wave1_adj = np.concatenate([[0], wave1, [0]])
t_wave1_adj = np.concatenate([[t_wave1[0]], t_wave1, [max(t_wave1)]])

wave2_adj = np.concatenate([[0], wave2, [0]])
t_wave2_adj = np.concatenate([[t_wave2[0]], t_wave2, [max(t_wave2)]])

# Create figure and axes
fig, axs = plt.subplots(1, 2, figsize=(12, 6), sharey=True)  # Two columns

colors = tc.tol_cset('muted')

# Plot original waveforms
axs[0].axhline(0, color='gray', linestyle='-', alpha=0.5)
axs[0].plot(t_wave1_adj * 1e6, -wave1_adj, label="Asymmetric & Bipolar", lw=1, color='k')

# now plot "Asymmetric monopolar", which is identical to the asymmetric biphasic waveform
axs[0].set_title("cathodic contact")
axs[0].set_xlabel("Time (μs)")
axs[0].set_ylabel("Unscaled amplitude")

# axs[0].legend()
axs[0].spines['top'].set_visible(False)
axs[0].spines['right'].set_visible(False)
# axs[0].legend(bbox_to_anchor=[3.4, 0.67])
# Plot polarity-flipped waveforms
axs[1].axhline(0, color='gray', linestyle='-', alpha=0.5, zorder=-1)

axs[1].plot(t_wave1_adj * 1e6, wave1_adj, label="Flipped Asymmetric Bipolar", lw=1, color='k')

axs[1].set_title("anodic contact")
axs[1].set_xlabel("Time (μs)")
axs[1].set_ylabel("")
# axs[1].legend()
axs[1].spines['top'].set_visible(False)
axs[1].spines['right'].set_visible(False)
plt.gcf().set_size_inches(4, 1.8)

# %%
