import matplotlib.pyplot as plt
import pandas as pd


inFiles = "final_spatial_trials_"

positions = []
lrs = []

#Load the files
for i in range(1, 45):
    with open(inFiles + str(i) + ".txt", "r") as file:
        for line in file:
            entries = line.split("\t")
            if entries[0] != "location":
                position = entries[0]
                lr = entries[1]
                positions.append(int(position))
                lrs.append(float(lr))
                
dataframe = pd.DataFrame()

dataframe['positions'] = positions
dataframe['lr'] = lrs

scatter = dataframe.scatter(x='positions', y='lr')
scatter.ticklabel_format(style='plain')
scatter.plot()
plt.show()
