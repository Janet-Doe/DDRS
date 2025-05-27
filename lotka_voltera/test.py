from math import *

def H(r,th,a,b,c,d,p0,q0):
    return a*(r*sin(th)-exp(r*sin(th)))+c*(r*cos(th)-exp(r*cos(th)))-a*(q0-exp(q0))-c*(p0-exp(p0))

def periode(a, b, c, d, x0, y0, N=1000):
    '''
    a, b, c, d : paramètres strictement positifs du modèle de Lotka-Voltera de base
    x0, y0 : conditions initiales
    N : nombre de points pour le calcul de l'inégrale, égal à 1000, par défaut
    '''
    if (x0==c/d) and (y0==a/b):
        raise("Ce calcul ne peut marcher que si on part de l'équilibre")
    p0=log(d*x0/c)
    q0=log(b*y0/a)
    g = array()
    g=zeros(1,N)
    