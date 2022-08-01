import matplotlib.pyplot as plt
import csv

file = open('submit/n_sims/0_0_2_0/data/outputs/Vm_original.dat', 'r')
text = file.readlines()
file.close()

ASCENT_x, ASCENT_y = [], []
for line in text:
    x_y = line.split(' ')
    ASCENT_x.append(float(x_y[0]))
    ASCENT_y.append(float(x_y[1]))

file = open('submit/n_sims/0_0_2_0/data/outputs/Vm_kna.dat', 'r')
text = file.readlines()
file.close()

kna_x, kna_y = [], []
for line in text:
    x_y = line.split(' ')
    kna_x.append(float(x_y[0]))
    kna_y.append(float(x_y[1]))

file = open('submit/n_sims/0_0_2_0/data/outputs/Vm_kdrTiger.dat', 'r')
text = file.readlines()
file.close()

kdrTiger_x, kdrTiger_y = [], []
for line in text:
    x_y = line.split(' ')
    kdrTiger_x.append(float(x_y[0]))
    kdrTiger_y.append(float(x_y[1]))

file = open('submit/n_sims/0_0_2_0/data/outputs/Vm_both.dat', 'r')
text = file.readlines()
file.close()

both_x, both_y = [], []
for line in text:
    x_y = line.split(' ')
    both_x.append(float(x_y[0]))
    both_y.append(float(x_y[1]))

file = open('submit/n_sims/0_0_2_0/data/outputs/Vm_both_Ra.dat', 'r')
text = file.readlines()
file.close()

both_Ra_x, both_Ra_y = [], []
for line in text:
    x_y = line.split(' ')
    both_Ra_x.append(float(x_y[0]))
    both_Ra_y.append(float(x_y[1]))

paper_x, paper_y = [], []
with open("validation/vm/literature/TIGERHOLM.csv", 'r') as file:
    csvreader = csv.reader(file)
    for row in csvreader:
        paper_x.append(float(row[0]))
        paper_y.append(float(row[1]))
paper_x = [i - min(paper_x) for i in paper_x]

peak_kdrTiger = kdrTiger_x[kdrTiger_y.index(max(kdrTiger_y))]
peak_paper = paper_x[paper_y.index(max(paper_y))]
paper_x = [i - abs(peak_kdrTiger-peak_paper) for i in paper_x]

zoom = 1
plt.style.use('seaborn-white')
fig = plt.figure()
ax1 = fig.add_subplot(111)

ax1.set_ylabel('v (mV)', fontsize=14)
ax1.set_xlabel('t (ms)', fontsize=14)
ax1.spines["top"].set_visible(False)
ax1.spines["right"].set_visible(False)
if zoom:
    ax1.plot(kdrTiger_x, kdrTiger_y, c='b', label='Adjusted kdrTiger', alpha=0.5)
    ax1.plot(both_x, both_y, c='r', label='Adjusted kdrTiger & Kna', alpha=0.5)
    ax1.plot(paper_x, paper_y, c='black', label='Tigerholm paper', alpha=0.5)
    ax1.plot(both_Ra_x, both_Ra_y, c='orange', label='Adjusted kdrTiger, Kna, & Ra', alpha=0.5)
    ax1.legend(loc='upper right')
    ax1.set_xlim(8, 18)
    ax1.set_ylim(-60, 45)
else:
    ax1.plot(ASCENT_x, ASCENT_y, c='b', label='Original ASCENT', alpha=0.5)
    ax1.plot(kna_x, kna_y, c='r', label='Adjusted Kna', alpha=0.5)
    ax1.plot(kdrTiger_x, kdrTiger_y, c='y', label='Adjusted kdrTiger', alpha=0.5)
    ax1.plot(both_x, both_y, c='g', label='Adjusted kdrTiger & Kna', alpha=0.5)
    ax1.plot(both_Ra_x, both_Ra_y, c='orange', label='Adjusted kdrTiger, Kna, & Ra', alpha=0.5)
    ax1.plot(paper_x, paper_y, c='black', label='Tigerholm paper', alpha=0.5)
    ax1.legend(loc='upper right')
    ax1.set_xlim(0, 50)
plt.show()
