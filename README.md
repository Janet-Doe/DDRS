---
title: "Rendu DDRS"
author: "Melvyn, Lilas, Simon"
output: 
  pdf_document:
    latex_engine: pdflatex
header-includes:
  - \usepackage{amsmath}
---

# DDRS
BAUVENT Melvyn, GRENIER Lilas, PRIBYLSKI Simon

## Modèle de Malthus, modèle de Verhulst

### Modèle de Malthus 
> Voir le code [ici](https://github.com/Janet-Doe/DDRS/tree/main/malthus).  

Afin d'afficher une analyse par le modèle de Malthus d'une population selon le temps, il est nécessaire de définir deux arrays correspondant respectivement aux années et à la population à ces années, et grâce aux librairies `numpy` et `matplotlib.pyplot` il suffit de faire : 

```py
coef = np.polyfit(years, np.log(population), 1)
poly1d_fn = np.poly1d(coef)
fig, ax = plt.subplots()
ax.plot(years, np.log(population), marker='+', linestyle='')
ax.plot(years, poly1d_fn(years), '-k')
ax.set_title('Population par année')
ax.set_xlabel('Année')
ax.set_ylabel('Population (M)')
plt.show()
```

### Modèle de Verhulst

> Voir le code [ici](https://github.com/Janet-Doe/DDRS/tree/main/verhulst).  

Consultez le fichier `verhulst57.py` pour un exemple d'utilisation de la régression linéaire.


## Modèle de Lotka-Voltera
> Voir le code [ici](https://github.com/Janet-Doe/DDRS/tree/main/lotka_voltera).

Pour afficher les courbes selon le modèle Lotka-Voltera, complétez le fichier main du dossier associé par des appels à la fonction `lotka_volterra_bis` avec les valeurs de a, b, c, d, x0 et y0 choisies.
