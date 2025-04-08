from logistic_regression import *
import numpy as np
import matplotlib.pyplot as plt

# Data
years=(1790, 1800, 1810, 1820, 1830, 1840, 1850, 1860, 1870, 1880, 1890, 1900, 1910, 1920, 1930, 1940, 1950)
population=(3.929, 5.308, 7.240, 9.638, 12.866, 17.069, 23.192, 31.443, 38.558, 50.156, 62.948, 75.995, 91.972, 105.711, 122.775, 131.669, 150.697)

t = np.array(years)
y = np.array(population)
anfin = None
P=1000

# Modelisation
K, y0, a, t0, resnorm, resnormc, yfin = regression_logistique(t, y)
popt = K, y0, a

coef = np.polyfit(years, np.log(population), 1)
poly1d_fn = np.poly1d(coef)


# Plot
fig, ax = plt.subplots()
ax.plot(years, population, marker='o', linestyle='-')
tc = np.linspace(min(t), max(t) if anfin is None else anfin, P)
plt.plot(tc, logistic_function(tc, *popt, t0), label='Modèle de Verhulst')
ax.plot(years, np.exp(poly1d_fn(years)), label='Modèle de Malthus')
ax.set_title('Population of USA (1790-1950)')
ax.set_xlabel('Year')
ax.set_ylabel('Population (millions)')
plt.show()

# Affichage des résultats
print(f"K estimé: {K}")
print(f"y0 estimé: {y0}")
print(f"a estimé: {a}")
print(f"t0 estimé: {t0}")
