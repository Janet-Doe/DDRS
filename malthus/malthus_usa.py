import numpy as np
import matplotlib.pyplot as plt

# Data
years=(1790, 1800, 1810, 1820, 1830, 1840, 1850, 1860, 1870, 1880, 1890, 1900, 1910, 1920, 1930, 1940, 1950)
population=(3.929, 5.308, 7.240, 9.638, 12.866, 17.069, 23.192, 31.443, 38.558, 50.156, 62.948, 75.995, 91.972, 105.711, 122.775, 131.669, 150.697)

# Modelisation
coef = np.polyfit(years, np.log(population), 1)
poly1d_fn = np.poly1d(coef)

# Plot
fig, ax = plt.subplots()
ax.plot(years, population, marker='o', linestyle='-')
ax.plot(years, np.exp(poly1d_fn(years)), '--k')
ax.set_title('Population of Canada (1851-2006)')
ax.set_xlabel('Year')
ax.set_ylabel('Population (millions)')
plt.show()