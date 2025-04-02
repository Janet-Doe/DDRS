function trace_droite_fenetre(coeff,r)

% tracé d'une droite dans la fenêtre d'un graphique déjà crée.
%
% *********************************************************
%  trace_droite_fenetre(coeff,r) :  
%  trace l'intersection de la droite d'équation aX+bY+C=0 où coeff=[a,b,c]
%  avec le pavé de coordonnées [X,Y] qui correspond à la fenêtre d'un
%  graphique déjà crée.
%  ATTENTION : les instructions graphiques appelant cette fonction doivent
%  être encadré par "hold on" et "hold off"
%  * Variables d'entrées : 
%     *  coeff :  tableau de type (1,3) représentant les trois coefficients
%     de l'équation de la droite.
%  * Variables d'entrée optionnelle :
%     * r  chaîne de caractère pour l'option du plot
%
% ************ Fonctions auxiliaires utilisées ************
%
%    aucune
%
% *********************************************************
%
% 2008 by          Jérôme BASTIEN
%                  Université Claude Bernard Lyon I, UFRSTAPS, Laboratoire CRIS, Villeurbanne
%                  E-Mail : jerome.bastien@univ-lyon1.fr
%

nin=nargin;
h=gca;
X=get(h,'XLim');
Y=get(h,'YLim');
xm=X(1);
xM=X(2);
ym=Y(1);
yM=Y(2);
if max(abs(coeff(1:2)))==0
    error('les deux coefficients de la droite sont nuls');
end
a=coeff(1);
b=coeff(2);
c=coeff(3);
test=1;
U=[];
V=[];
if b==0
    test=0;
    x=-c/a;
    if ((xm<=x)&(x<=xM))
        U=x([1 1]);
        V=Y;
    end
end
if test&(a==0)
    test=0;
    y=-c/b;
    if ((ym<=y)&(y<=yM))
        U=X;
        V=y([1 1]);
    end
end
if test
    alphaa=-a/b;
    betaa=-c/b;
    y=alphaa*X+betaa;
    if y(1)>yM
        if y(2)<=yM
            if y(2)<ym
                x=(Y-betaa)/alphaa;
                U=x;
                V=Y;
            else
                x=(yM-betaa)/alphaa;
                U=[x,xM];
                V=[yM,y(2)];
            end
        end
    else
        if y(1)<ym
            if y(2)>yM
                x=(Y-betaa)/alphaa;
                U=x;
                V=Y;
            else
                if y(2)>=ym
                    x=(ym-betaa)/alphaa;
                    U=[x,xM];
                    V=[ym,y(2)];
                end
            end
        else
            if y(2)>yM
                x=(yM-betaa)/alphaa;
                U=[xm,x];
                V=[y(1),yM];
            else
                if y(2)<ym
                    x=(ym-betaa)/alphaa;
                    U=[xm,x];
                    V=[y(1),ym];
                else
                    U=X;
                    V=y;
                end
            end
        end
    end
end
if nin==1
    plot(U,V);
else
    plot(U,V,r);
end