import glob
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit


def gaussian(x, a, x0, sigma):
    return a*np.exp(-(x-x0)**2/(2*sigma**2))

file_list = glob.glob('*.txt')
location_all = []
LR_all = []
alpha_all = []
means = []
variances = []
amplitudes = []
popts = []

for file_name in file_list:
    print(file_name)
    df = pd.read_csv(file_name, sep='\t', header=0)
    location = df['location']
    LR = df['LR']
    alpha = df['alpha']
    location_all.append(location)
    LR_all.append(LR)
    alpha_all.append(alpha)

    # Fit a Gaussian distribution to the LR data
    popt, _ = curve_fit(gaussian, location, LR,maxfev=10000)
    popts.append(popt)
    amplitudes.append(popt[0])
    means.append(popt[1])
    variances.append(popt[2]**2)

# LR_all vs. position_all scatter plot with Gaussian fit
fig, (ax1) = plt.subplots(nrows=1, ncols=1, sharex=True)

for location, LR in zip(location_all, LR_all):
    ax1.scatter(location, LR, alpha=0.5)

ax1.plot(location_all[0], gaussian(location_all[0], *popt), 'r-')
ax1.set_ylabel('LR')
ax1.set_title('well_mixed LR vs. Position')

plt.show()

# Print means, amplitudes, and variances
print('Mean, amplitude, and variance of Gaussian fit to LR data:')
print(pd.DataFrame({'Mean': means, 'Amplitude': amplitudes, 'Variance': variances}))

# Calculate mean of means and variance of the mean
mean_of_means = np.mean(means)
variance_of_means = np.sum((means - mean_of_means)**2)/len(file_list)
mean_of_amplitudes = np.mean([popt[0] for popt in popts])


print(f'Mean of means: {mean_of_means}')
print(f'Mean of amplitudes: {mean_of_amplitudes}')
print(f'Variance of the mean: {variance_of_means}')
