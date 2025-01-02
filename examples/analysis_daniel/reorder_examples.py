"""Created on Mon Aug 28 16:02:12 2023.

@author: dpm42
"""

import math


def compute_reorder_cost(order1, order2):
    assert set(order1) == set(order2)
    assert len(set(order1)) == len(order1)
    assert len(set(order2)) == len(order2)
    sumdist = 0
    for item in order1:
        distance = abs(order1.index(item) - order2.index(item))
        sumdist += distance
    if len(order1) % 2 == 0:
        maxlen = sum(range(len(order1))[int(math.ceil(0.5 * len(order1))) :]) * 2 / len(order1)
    else:
        maxlen = (
            range(len(order1))[int(math.floor(0.5 * len(order1)))]
            + (sum(range(len(order1))[int(math.ceil(0.5 * len(order1))) :])) * 2
        ) / len(order1)
    return sumdist / len(order1) / maxlen


# compute maximal reorder cost
input1 = list(range(2))
input2 = list(reversed(range(2)))
print(compute_reorder_cost(input1, input2))
# plot change in the above as the length of the list increases
import matplotlib.pyplot as plt
import numpy as np

x = np.arange(2, 100, 5)
y = [compute_reorder_cost(list(range(i)), list(reversed(range(i)))) for i in x]
plt.scatter(x, y)
plt.title('Reorder cost for reversed lists')
plt.xlabel('Length of list')
plt.ylabel('AR')
plt.show()
# compute reorder cost for a random list
import random

input1 = list(range(100))
# rearrange the list without replacement
input2 = random.sample(input1, len(input1))
print(compute_reorder_cost(input1, input2))
# plot change in the above as the length of the list increases
# do 20 iterations for each length and plot the mean and standard deviation
import seaborn as sns

x = np.arange(2, 100, 5)
y = []
for i in x:
    y.append([compute_reorder_cost(list(range(i)), random.sample(list(range(i)), i)) for _ in range(20)])
y = np.array(y)
plt.errorbar(x, y.mean(axis=1), yerr=y.std(axis=1), fmt='o')
plt.title('Reorder cost for random lists')
plt.xlabel('Length of list')
plt.ylabel('AR')
plt.show()

# now do the same if the list is "shifted" by half its length
input2 = input1[int(len(input1) / 2) :] + input1[: int(len(input1) / 2)]
print(compute_reorder_cost(input1, input2))
# plot change in the above as the length of the list increases

x = np.arange(2, 100, 5)
y = []
for i in x:
    y.append(
        [compute_reorder_cost(list(range(i)), list(range(int(i / 2), i)) + list(range(int(i / 2)))) for _ in range(20)]
    )
y = np.array(y)
plt.errorbar(x, y.mean(axis=1), yerr=y.std(axis=1), fmt='o')
plt.title('Reorder cost for shifted lists')
plt.xlabel('Length of list')
plt.ylabel('AR')
plt.show()
