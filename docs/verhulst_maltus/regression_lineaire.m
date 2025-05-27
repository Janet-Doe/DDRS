function [a,b,r]=regression_lineaire(x,y)

% regression_lineaire : calcul de droite de moindres carrées.
%
% x et y sont deux vecteurs de même taille.
% [a,b,r]=regression_lineaire(x,y) renvoie a : le pente de la droite,
%                                          b : l'ordonnée à l'origine,
%                                          r : la correlation
%
% 
% 
% ************ Fonctions auxiliaires utilis�es ************
%
%       aucune
%
% *********************************************************
%

r=corrcoef(x,y);
r=r(1,2);
[p,s]=polyfit(x,y,1);
a=p(1);
b=p(2);



