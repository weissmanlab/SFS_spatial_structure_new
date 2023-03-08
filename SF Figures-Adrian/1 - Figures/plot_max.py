import matplotlib.pyplot as plt
import pandas as pd


inFiles = "final_spatial_trials_"

maximumPos = []
#Load the files
for i in range(1, 45):
    currMax = -10000
    currPos = -10
    with open(inFiles + str(i) + ".txt", "r") as file:
        for line in file:
            entries = line.split("\t")
            if entries[0] != "location":
                position = int(float(entries[0]))
                lr = float(entries[1])
                if lr > currMax:
                    currMax = lr
                    currPos = position
        maximumPos.append(currPos)
                
dataframe = pd.DataFrame()

dataframe['positions'] = maximumPos
hist = dataframe.hist()
plt.show()
