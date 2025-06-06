# DDRS
BAUVENT Melvyn, GRENIER Lilas, PRIBYLSKI Simon

## Modèle de Malthus, modèle de Verhulst
Modèle de Malthus : [ici](malthus).  
Modèle de Verhulst : [ici](verhulst).  

### Exercice 1.1

#### 1a

La solution de l'équation différentielle est :

$$
    \begin{cases}
      N(t)=N_0 e^{r(t-t_0)}\\
    \end{cases}
$$

On a bien :

$$
    \begin{cases}
      N'(t)=N_0 r e^{r(t-t_0)} = r N(t) \\
      N(t_0) = N_0 e^0 = N_0
    \end{cases}
$$

#### 1b

$$
    \text{Pour r > 0 :} \lim_{x \to +\infty} N(t) = +\infty\\
    \text{Pour r < 0 :} \lim_{x \to +\infty} N(t) = 0
$$

$$
    N(t) > 0\\
    \text{Pour r > 0 : } N'(t) = rN(t) > 0\\
    \text{Pour r < 0 : } N'(t) = rN(t) < 0
$$

#### 1c

$$
    \tau(t)=\frac{N(t+h)-N(t)}{h}\frac{1}{N(t)}\\
    \lim_{h \to 0} \tau(t) = \frac{N'(t)}{N(t)} = r
$$

#### 2a

Avec b = 0, on a :

$$
    N'(t)=aN(t)
$$

On retrouve le modèle de Malthus

#### 2b(i)

En prenant

$$
    N(t) = K\\
    N(t_0)=K
$$

On obtient :

$$
    \begin{cases}
        N'(t)=0\\
        a N(t) (1 - \frac{N(t)}{K}) = a K (1 - 1) = 0
    \end{cases}
$$

#### 2b(ii)

$$
    \text{On suppose que N(t) > 0, on pose v =} \frac{1}{N}\\
    v' + av = (\frac{1}{N})'+ \frac{1}{N} = \frac{-N}{N^2} = \frac{-aN(1-\frac{N}{K})}{N^2} + \frac{a}{N} = \\
    \frac{-a(1-\frac{N}{K})}{N} + \frac{a}{N} = \frac{a(\frac{N}{K})}{N} = \frac{a}{K}
$$

#### 2b(iii)

Solution de l'équation homogène :

$$
v_h(t) = Ce ^ {-a(t-t_0)}
$$

Solution particulière :

$$
v_s(t) = \frac{1}{K}
$$







## Modèle de Lotka-Voltera
Voir le code [ici](lotka_voltera).

### Point d'équilibre 
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

### Variations des paramètres initiaux.
x_0

Paramètres :   
$a$ : nombre de prédateurs dans la population.  
$c$ : taux de croissance des proies.   
$b$ : $d$x de croissance des prédateurs.   
$x_0$ : nombre de proies au moment $t_0$.  
$y_0$ : nombre de prédateurs au moment $t_0$.  

Résultats :  
$x$ : nombre de proies dans la population.  
$y$ : nombre de prédateurs dans la population.   

Dans l'exemple du cours, ces paramètres prenaient les valeurs : 

$$ 
    a = 0, 734 545 7 \\
    b = d = 2, 583 453 056 × 10^{−5} \\
    c = 0, 598 414 3 \\
    x_0 = 79 549, 777 045 0 \\
    y_0 = 35 216, 498 070 0
$$

Les graphes obtenus à partir de ces valeurs sont les suivants :
`b` = `d` :  efficacité de la prédation. 
``
#### Variation de a : 
e la même, mais la moyenne de population de prédateurs croît.  
Voir schéma, pour les valeurs $ a = [1, 2, 3, 4, 5] $ :
`````b
#### Variation de b et d : 
lace s

Voir schéma, pour les valeurs $b = d = [1, 1.5, 2, 2.5, 3] * 10^{-5} $ :
````c
#### Variation de c : 

Voir schéma, pour les valeurs $c = [1, 2, 3, 4, 5]$ :
``````
#### Variation simultanée de a et c : 

On peut également faire varier les paramètres a et c simultanément.

#### Variation de $x_0$ et $y_0$ : 
iales, on obtient un schéma de forme similaire, de même point d'équilibre (puisque $x_0$ et $y_0$ n'entrent pas en compte dans le calcul du point d'équilibre tant qu'ils sont >0) mais d'envergure différente.  
Schéma pour valeurs $x_0 = y_0 = [1, 2, 3, 4, 5]$ :

#### Conservation du point d'équilibre


/todoOn cherche à garder le même point d'équilibre $(\frac{c}{d}, \frac{a}{b})$ en changeant les paramètres : 
