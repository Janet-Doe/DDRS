import numpy as np
import matplotlib.pyplot as plt

# Data
years = (1851, 1861, 1871, 1881, 1891, 1901, 1911, 1921, 1931, 1941, 1951, 1956, 1961, 1966, 1971, 1976, 1981, 1986, 1991, 1996, 2001, 2006)
population = (2.436297, 3.229633, 3.737257, 4.381256, 4.932206, 5.418663, 7.221662, 8.800429, 10.376379, 11.506655, 14.009429, 16.080791, 18.238247, 20.014880, 21.568305, 22.992595, 24.343177, 25.309330, 27.296856, 28.846758, 30.007094, 31.612897)

# Modelisation
coef = np.polyfit(years, population, 1)
poly1d_fn = np.poly1d(coef) 

# Plot
fig, ax = plt.subplots()
ax.plot(years, population, marker='o', linestyle='-')
ax.plot(years, poly1d_fn(years), '--k')
ax.set_title('Population of Canada (1851-2006)')
ax.set_xlabel('Year')
ax.set_ylabel('Population (millions)')
plt.show()