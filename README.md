# DDRS
BAUVENT Melvyn, GRENIER Lilas, PRIBYLSKI Simon

## Modèle de Malthus 
Voir [ici](malthus).

## Modèle de Verhulst
Voir [ici](verhulst).

## Modèle de Lotka-Voltera
Voir le code [ici](lotka_voltera).

### Exercice 2.1

#### Partie 1
On veut montrer que :
$$
    \text{(1) : }
    x'(t)=y'(t)=0 \Leftrightarrow 
    \begin{cases}
      x(t)=0 \text{ ou } y(t)= \frac{a}{b}\\
      y(t)=0 \text{ ou } x(t)= \frac{c}{d}
    \end{cases}       
$$ 

Pour cela, on souhaite utiliser les équations (43) du support de cours, avec $a$, $b$, $c$, $d$, $x_0$ et $y_0 > 0$ :

$$
    \begin{cases}
        x'(t) = ax(t) - bx(t)y(t) \\
        y'(t) = -cy(t) + dx(t)y(t) \\
        x(0) = x_0 \\
        y(0) = y_0
    \end{cases}   
$$ 

En utilisant le côté gauche de l'expression (1), on obtient : 

$$
    \begin{cases}
        ax(t) = bx(t)y(t) \\
        cy(t) = dx(t)y(t) \\
        x(0) = x_0 \\
        y(0) = y_0
    \end{cases}   

\Leftrightarrow 

    \begin{cases}
        x(t) = 0 \text{ ou } y(t) = \frac{a}{b} \\
        y(t) = 0 \text{ ou } x(t) = \frac{c}{d} \\
        x(0) = x_0 \\
        y(0) = y_0
    \end{cases}   

$$ 

De plus, avec les conditions initiales $x_0$ et $y_0 > 0$, les conditions $
    \begin{cases}
        x(0) = x_0 \\
        y(0) = y_0
    \end{cases}   
$ nous permettent d'affirmer que nous ne pouvons avoir ni $x(t) = 0$ ni $y(t) = 0$. 

Ainsi, nous obtenons la relation :

$$
    \text{(2) : }
    x'(t)=y'(t)=0 \Leftrightarrow 
    \begin{cases}
        x(t) = \frac{c}{d} \\
        y(t) = \frac{a}{b} \\
        x(0) = x_0 \\
        y(0) = y_0
    \end{cases}  
$$

Le seul point d'équilibre possible est donc le point $(\frac{c}{d}, \frac{a}{b})$.

#### Partie 2

