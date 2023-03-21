# mc.py : Monte Carlo calculation on ising model with triangle lattice

import os
num_thread = 1
os.environ['OMP_NUM_THREADS'] = str(num_thread)
os.environ['OPENBLAS_NUM_THREADS'] = str(num_thread)

import re
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

with open('test.txt', 'r') as f: data = np.genfromtxt(f)

fig, ax = plt.subplots(1, 2, figsize=(10, 5), constrained_layout=True)
for i, title in zip(range(data.shape[1]), ['Energy', 'Magnetization']):
	ax[i].plot(data[:, 0], data[:, i+1], marker='.')
	ax[i].set_title(title)
	ax[i].grid()
fig.supxlabel('Temperature(T)')
fig.savefig('diagram/mag.png')
plt.show()
