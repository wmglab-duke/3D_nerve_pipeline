import matplotlib.pyplot as plt
import pickle
import os

# file = open('extracellular_stim/vm/literature/', 'r')
# text = file.readlines()
# file.close()
# 
# HOC_x, HOC_y = [], []
# for line in text:
#     x_y = line.split(' ')
#     HOC_x.append(float(x_y[0]))
#     HOC_y.append(float(x_y[1]))
    
with open('extracellular_stim/vm/Python/Vm_TIGERHOLM', 'rb') as f:
    Python_vals = pickle.load(f)
old_x, old_y = Python_vals[0], Python_vals[1]

with open('../Vm_TIGERHOLM', 'rb') as f:
    Python_vals = pickle.load(f)
new_x, new_y = Python_vals[0], Python_vals[1]

plt.style.use('seaborn-white')
fig = plt.figure()
ax1 = fig.add_subplot(111)

ax1.set_ylabel('v (mV)', fontsize=14)
ax1.set_xlabel('t (ms)', fontsize=14)
ax1.spines["top"].set_visible(False)
ax1.spines["right"].set_visible(False)
ax1.plot(new_x, new_y, c='orange', label='New', alpha=0.5)
ax1.plot(old_x, old_y, c='blue', label='Old', alpha=0.5)
ax1.legend(loc='upper right')

plt.show()