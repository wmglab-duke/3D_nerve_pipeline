import numpy as np
from matplotlib import pyplot as plt

fascs = [4.57,4.6,11.91,8.97,12.15,6.36]
branch = [68,104,274,136,195,128]
plt.scatter(fascs,branch)
plt.title('Fascicle count vs branch count')
plt.xlabel('Fascicle count')
plt.ylabel('Branch count')
#calculate correlation and add to plot
correlation = np.corrcoef(fascs,branch)
#round
correlation = np.round(correlation, decimals=2)
plt.text(0.2, 0.8, f'Correlation: {correlation[0,1]}', horizontalalignment='center', verticalalignment='center', transform=plt.gca().transAxes)
#add regression line
m, b = np.polyfit(fascs, branch, 1)
plt.plot(fascs, m*np.array(fascs) + b,'r')
plt.show()