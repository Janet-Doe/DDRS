# Copyright 2025 BAUVENT Melvyn, GRENIER Lilas, PRIBYLSKI Simon

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit


def logistic_function(t, K, y0, a, t0):
    return K / (1 + (K / y0 - 1) * np.exp(-a * (t - t0)))


def regression_logistique(t, y, modmal=1, P=1000, yc=None, ylabelp='', MaxFunEvals=None, cho=1, iter=0,
                          anfin=None, tansym=0):
    if np.any(y <= 0):
        raise ValueError("Il existe des valeurs négatives ou nulles de y")

    # Initialisation avec la méthode de Perrin
    X = y[:-1]
    Y = (y[1:] / X - 1) / np.diff(t)
    p = np.polyfit(X, Y, 1)
    a0 = p[1]
    K0 = -a0 / p[0]
    y00 = y[0]

    if a0 <= 0 or K0 <= 0:
        raise ValueError("a0 <= 0 ou K0 <= 0")

    t0 = t[0]
    x0 = [K0, y00, a0]

    if a0 <= 0 or K0 <= 0 or y00 <= 0 or y00 >= K0:
        raise ValueError("a0 <= 0 ou K0 <= 0 ou y00 <= 0 ou y00 >= K0")

    # Ajustement avec curve_fit
    maxfev = MaxFunEvals if MaxFunEvals is not None else 10000  # Correction de NoneType
    try:
        popt, pcov = curve_fit(lambda t, K, y0, a: logistic_function(t, K, y0, a, t0), t, y, p0=x0, maxfev=maxfev)
    except RuntimeError:
        popt, pcov = curve_fit(lambda t, K, y0, a: logistic_function(t, K, y0, a, t0), t, y, p0=x0, bounds=(0, np.inf),
                               maxfev=maxfev)

    K, y0, a = popt

    if a <= 0 or K <= 0 or y0 <= 0:
        raise ValueError("a <= 0 ou K <= 0 ou y0 <= 0")

    resnorm = np.sum((y - logistic_function(t, *popt, t0)) ** 2)
    resnormc = None
    yfin = logistic_function(anfin, *popt, t0) if anfin is not None else None

    return K, y0, a, t0, resnorm, resnormc, yfin


# Code de test
if __name__ == "__main__":
    # Données réelles de la population des USA
    t = np.array([1790, 1800, 1810, 1820, 1830, 1840, 1850, 1860, 1870, 1880, 1890, 1900, 1910, 1920, 1930, 1940, 1950])
    y = np.array(
        [3.929, 5.308, 7.240, 9.638, 12.866, 17.069, 23.192, 31.443, 38.558, 50.156, 62.948, 75.995, 91.972, 105.711,
         122.775, 131.669, 150.697])

    # Exécution de la régression logistique
    K, y0, a, t0, resnorm, resnormc, yfin = regression_logistique(t, y)

    # Affichage des résultats
    print(f"K estimé: {K}")
    print(f"y0 estimé: {y0}")
    print(f"a estimé: {a}")
    print(f"t0 estimé: {t0}")