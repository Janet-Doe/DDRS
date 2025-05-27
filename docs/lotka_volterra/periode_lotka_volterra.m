function [T,er]=periode_lotka_volterra(a,b,c,d,x0,y0,N)

% Détermination de la période du modèle de base de Lotka_Voltera
%
% *********************************************************
% [T,er]=periode_lotka_volterra(a,b,c,d,x0,y0,N) :
%
% Pour le modèle de base
%        x'(t) =x*(a*x-b*y)
%        y'(t) =y*(-c*y+d*x)
%	* Entrées :
%     * a,b,c,d : paramètres strictement positifs du modèle de
%       Lotka-Voltera de base
%     * x0,y0 : conditions initiales
%	* Entrée optionnelle :
%     * N : nombre de points pour le calcul de l'inégrale, égal à 1000, par
%       défaut
%	* Sortie :
%     * T : la période
%	* Sortie facultative :
%     * er : indicateur d'erreur pour la résolution des équations non
%     linéaires.
%
% ************ Fonctions auxiliaires utilisées ************
%
%    aucune
%
% *********************************************************
%
% (c) 2025 by      Jérôme BASTIEN,
%                  LIBM, Polytech, Université Lyon 1
%                  E-Mail : jerome.bastien@univ-lyon1.fr
%
% Exemples
%    T=periode_lotka_volterra(1,1,1,1,0.5,0.5);
%    T=periode_lotka_volterra(3,1,2,1,1,2); %celui de G. Vial

% Ce calcul ne peut marcher si on part de l'équilibre !!!
if x0==c/d&&y0==a/b
    error('Ce calcul ne peut marcher si on part de l''équilibre');
end
if nargin<=6||~isempty(N)
    N=1000;
end
p0=log(d*x0/c);
q0=log(b*y0/a);
H=@(r,th,a,b,c,d,p0,q0) a*(r*sin(th)-exp(r*sin(th)))+c*(r*cos(th)-exp(r*cos(th)))-a*(q0-exp(q0))-c*(p0-exp(p0));
g=zeros(1,N);
if nargout>=2
    er=g;
end
[th0,r0] =cart2pol(p0,q0);
auxi=r0;
g(1)=r0;
th=th0+linspace(0,2*pi,N);
if nargout>=2
    er(1)=abs(H(r0,th0,a,b,c,d,p0,q0));
end
for k=2:N
    auxib=fzero(@(r) H(r,th(k),a,b,c,d,p0,q0),auxi);
    if nargout>=2
        er(k)=abs(H(auxib,th(k),a,b,c,d,p0,q0));
    end
    g(k)=auxib;
    auxi=auxib;
end
h=(2*pi)/(N-1);
sinth=sin(th);
costh=cos(th);
T=-h*trapz(g./(a*sinth.*(1-exp(g.*sinth))+c*costh.*(1-exp(g.*costh))));
if nargout>=2
    er=max(abs(er));
end