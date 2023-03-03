import glob
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

file_list = glob.glob('*.txt')

counts_all = []
positions_all = []

for file_name in file_list:
    print(file_name)
    print(file_name)
    df = pd.read_csv(file_name, sep='\t', header=0)
    positions = df['position']
    counts = df['x']
    positions_all.append(positions)
    counts_all.append(counts)

fig, ax = plt.subplots()
positions_all = np.concatenate(positions_all)
counts_all = np.concatenate(counts_all)
ax.hist(positions_all, bins=100, alpha=0.5, weights=(counts_all/(1/100)/(15*100))) #weights=counts/window size / simple size

ax.set_xlabel('Position')
ax.set_ylabel(' polymorphic sites dentiy ')
ax.set_title('Distribution of Allele Frequencies in a 10x10 Structured Population')

plt.show()
