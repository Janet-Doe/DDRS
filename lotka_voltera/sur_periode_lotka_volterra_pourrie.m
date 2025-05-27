function [x,y,t,er]=sur_periode_lotka_volterra_pourrie(a,b,c,d,x0,y0,T,h)

% D�termination des solutions du mod�le de base de Lotka_Voltera sur une p�riode
%
% *********************************************************
% [X,Y]=sur_periode_lotka_volterra_pourrie(a,b,c,d,x0,y0,T,h) :
%
% Pour le mod�le de base
%        x'(t) =x*(a*x-b*y)
%        y'(t) =y*(-c*y+d*x)
%	* Entr�es :
%     * a,b,c,d : param�tres strictement positifs du mod�le de
%       Lotka-Voltera de base
%     * x0,y0 : conditions initiales
%     * T : p�riode d�j� calcul�e (par periode_lotka_volterra)
%	* Entr�e optionnelle :
%     * pas de discr�tisation, �gal � 1e-3 par d�faut.
%	* Sortie :
%     * x,y,t : valeurs de x et de y (et du temps t)
%	* Sortie optionnelle
%     * er : indicateur d'erreur pour la r�solution des �quations non
%     lin�aires.
%
%
% ************ Fonctions auxiliaires utilis�es ************
%
%    aucune
%
% *********************************************************
%
% (c) 2025 by      J�r�me BASTIEN,
%                  LIBM, Polytech, Universit� Lyon 1
%                  E-Mail : jerome.bastien@univ-lyon1.fr
%
% Exemple
%   a=1;b=1;c=1;d=1;x0=0.5;y0=0.5;
%   Tp=periode_lotka_volterra(a,b,c,d,x0,y0);
%   h=1e-3;
%   N=ceil(Tp/h);
%   h0=Tp/N;
%   [x,y,t,erb]=sur_periode_lotka_volterra_pourrie(a,b,c,d,x0,y0,Tp,h0);
%   plot(t,x,t,y);
%   figure;
%   [T,Xc,Yc,er]=lotka_volterra_complet(a,b,c,d,x0,y0,0,Tp,[],[],[],[],[],[],h0);% %% Attention, on a 
%   max(abs(t-T)) 
%   max([max(abs(Xc-x)) max(abs(Yc-y))]) 
%   max([max(abs(Xc-x)) max(abs(Yc-y))]) 
% %% ATTENTION, gard�e pour m�moire, cette fonction est fausse !!!
% %% car theta(t) v�rifie l'�quation diff�rentielle 
% %% theta'(t)=Frond(theta(t))
% %% theta(0)=theta_0
% %% et donc on n'a pas theta(t) fonction lin�aire de t comme je le croyais
% !!!
% N�anmoins, les lignes 3 � 6 du script donn�s en exemple  sont int�ressantes !!!
% On peut se contenter de calculer et de r�soudre sur une seule p�riode
% !!!!

if nargin<=7||isempty(h)
    h=1e-3;
end
% lignes adapt�e de periode_lotka_volterra
p0=log(d*x0/c);
q0=log(b*y0/a);
H=@(r,th,a,b,c,d,p0,q0) a*(r*sin(th)-exp(r*sin(th)))+c*(r*cos(th)-exp(r*cos(th)))-a*(q0-exp(q0))-c*(p0-exp(p0));
[th0,r0] =cart2pol(p0,q0);
t=(0:h:T).';
th=th0+2*pi*t/T;
N=length(t);
g=zeros(N,1);
if nargout>=3
    er=g;
end
auxi=r0;
g(1)=r0;
if nargout>=4
    er(1)=abs(H(r0,th0,a,b,c,d,p0,q0));
end
for k=2:N-1
    auxib=fzero(@(r) H(r,th(k),a,b,c,d,p0,q0),auxi);
    if nargout>=4
        er(k)=abs(H(auxib,th(k),a,b,c,d,p0,q0));
    end
    g(k)=auxib;
    auxi=auxib;
end
if nargout>=4
    er=max(abs(er));
end
g(end)=r0;
x=(c/d)*exp(g.*cos(th));
y=(a/b)*exp(g.*sin(th));