import glob
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
from scipy.stats import norm

# def gaussian(x, a, x0, sigma, b):
#     return -a * np.exp(-((x - x0) ** 2) / (2 * sigma ** 2)) + b

file_list = glob.glob('*.txt')
location_all = []
alpha_all = []
means = []
variances = []

for file_name in file_list:
    print(file_name)
    df = pd.read_csv(file_name, sep='\t', header=0)
    location = df['location']
    alpha = df['alpha']
    location_all.append(location)
    alpha_all.append(alpha)

    # Fit a negative Gaussian distribution to the alpha data
    # popt, _ = curve_fit(gaussian, location, alpha, p0=[1, 0.8, 1, 1], maxfev=10000)
    # means.append(popt[1])
    # variances.append(popt[2] ** 2)

# alpha_all vs. position_all scatter plot with negative Gaussian fit
fig, ax1 = plt.subplots(nrows=1, ncols=1, sharex=True)

for location, alpha in zip(location_all, alpha_all):
    ax1.scatter(location, alpha, alpha=0.5)

x = np.linspace(np.min(location_all), np.max(location_all), 100)
# ax1.plot(x, gaussian(x, *popt), 'r-', label='Negative Gaussian Fit')
ax1.set_ylabel('alpha')
ax1.set_title('well_mixed alpha vs. Position')

plt.show()

# # Print means and variances
# print('Mean and variance of negative Gaussian fit to alpha data:')
# print(pd.DataFrame({'Mean': means, 'Variance': variances}))

# # Calculate mean of means and variance of the mean
# mean_of_means = np.mean(means)
# variance_of_means = np.sum((means - mean_of_means) ** 2) / len(file_list)

# print(f'Mean of means: {mean_of_means}')
# print(f'Variance of the mean: {variance_of_means}')
