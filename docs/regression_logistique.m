function [K,y0,a,t0,resnorm,resnormc,yfin]=regression_logistique(t,y,modmal,P,yc,ylabelp,MaxFunEvals,cho,iter,grap,anfin,tansym)

% param�tres de K,y0,a,t0 de la regression logistique (mod�le de Verhulst) pour des donn�es (t,y)
%
% *********************************************************
%   [K,y0,a,t0,resnorm,resnormc,yfin]=regression_logistique(t,y,modmal,P,yc,ylabelp,MaxFunEvals,cho,iter,grap,anfin,tansym)
%
%	* Entr�es :
%     * t,y : deux tableaux ligne repr�sentant les donn�es ti,yi (dates,
%       population)
%	* Entr�es optionnelles :
%     * modmal entier dans {0,1} �gal � 1 par d�faut
%      * si 0 : pas d'affichage du mod�le de Malthus
%      * si 1 : affichage du mod�le de Malthus
%     * P : nombre de points pour le graphique continu, �gal � 1000 par
%     d�faut
%     * yc : autres mesures faites (avec les m�mes ordonn�es t, qui donnera
%       l'erreur en sortie resnormc), vide par d�faut.
%     * labelp : coefficient multiplicateur pour la l�gende "population" en
%     y, vide par d�faut.
%     * MaxFunEvals: vide par d�faut, pour l'option de lsqcurvefit
%     * cho : voir fonction fonction_logistique
%     * iter : entier dans {0,1}, �gal � 0 par d�faut.
%        affiche les r�sultats partiels de l'optimisation si iter=1.
%     * grap : entier dans {0,1}, �gal � 1 par d�faut.
%        produit le graphique si grap=1. Sur ce graphique, apparaissent :
%          * les donn�es mesur�es (ti,yi)
%          * la courbe correpondant au mod�le de Verhulst (voir la fonction
%            fonction_logistique)
%          * �ventuellement, la courbe correpondant au mod�le de Malthus
%          * la courbe correpondant au mod�le de Verhulst  dont les valeurs
%             sont dans yc
%     * anfin : ann�e de fin de graphique, qui n'est pas n�cessairement la
%       donn�e maximale de t. vide par d�faut. Si non vide, la valeur de la
%       population pour cette ann�e est renvoy�e dans yfin
%     * tansym : entier dans {0,1}, �gal � 0 par d�faut.
%        si =1, affiche la tangente au point d'inflexion.
%  * Sorties
%     * K,y0,a,t0 tel que yi=K/(1+(K/y0-1)*exp(-a*(t-t0)) au sens des
%       moindres carr�s, correspondant � la fonction logistique (voir
%       fonction_logistique)
%     Attention, t0 est arbitraire et sera choisi �gal � la premi�re valeur de t
%  * Sorties optionnelles
%     * resnorm : sommes des carr�s des �carts entre les donn�es et le
%       mod�le.
%     * resnormc : voir yc
%     * yfin : voir anfin
%
%
% ************ Fonctions auxiliaires utilis�es ************
%
%    fonction_logistique
%    regression_lineaire
%    trace_droite_fenetre
%
% *********************************************************
%
% (c) 2024 by      J�r�me BASTIEN,
%                  LIBM, Polytech, Universit� Lyon 1
%                  E-Mail : jerome.bastien@univ-lyon1.fr
%
% Exemples
% a)
%  donnee_usa;
%  [K,y0,a,t0]=regression_logistique(t,y,[],[],[],'M');
% b) (JB) voir fichiers elephant.matex et usa_verhulst.matex
%         dans \LINUX\enseignement\Polytech\2023_2024\DDRSP24
%

if nargin<=2||isempty(modmal)
    modmal=1;
end
if nargin<=3||isempty(P)
    P=1000;
end
if nargin<=4
    yc=[];
end
if nargin<=5
    ylabelp='';
end
if nargin<=6
    MaxFunEvals=[];
end
if nargin<=7||isempty(cho)
    cho=1;
end
if nargin<=8||isempty(iter)
    iter=0;
end
if nargin<=9||isempty(grap)
    grap=1;
end
if nargin<=10
    anfin=[];
end
if nargin<=11||isempty(tansym)
    tansym=0;
end
ind=y<=0;
if any(ind)
    error('Il existe des valeurs n�gatives ou nulles de y');
end
% Initialisation avec la M�thode de Perrin
X=y(1:end-1);
Y=(y(2:end)./X-1)./diff(t);
p=polyfit(X,Y,1);
a0=p(2);
K0=-a0/p(1);
y00=y(1);
if a0<=0||K0<=0
    error('a0<=0 ou K0<=0');
end
t0=t(1);
x0=[K0,y00,a0];
if a0<=0||K0<=0||y00<=0||y00>=K0
    error('a0<=0 ou K0<=0 ou y00<=0 ou y00>=K0');
end
% calcul avec lsqcurvefit
if iter
    if isempty(MaxFunEvals)
        options = optimset('Display','iter','TolFun',eps,'TolX',eps);
    else
        options = optimset('MaxFunEvals',MaxFunEvals,'Display','iter','TolFun',eps,'TolX',eps);
    end
else
    if isempty(MaxFunEvals)
        options =optimset('TolFun',eps,'TolX',eps);
    else
        options = optimset('MaxFunEvals',MaxFunEvals,'TolFun',eps,'TolX',eps);
    end
end
try
    lb=[];
    [x,resnorm,residual,flag]=lsqcurvefit(@(x,xdata) fonction_logistique(x,xdata,t0,cho),x0,t,y,lb,[],options);
catch
    lb=eps*ones(1,3);
    [x,resnorm,residual,flag]=lsqcurvefit(@(x,xdata) fonction_logistique(x,xdata,t0,cho),x0,t,y,lb,[],options);
end
disp('exitflag');
disp(flag);
K=x(1);
y0=x(2);
a=x(3);
if a<=0||K<=0||y0<=0
    error('a<=0 ou K<=0 ou y0<=0');
end
if ~isempty(yc)
    resnormc=sum((yc-y).^2);
    % On essaye de rappeler la fonction regression_logistique pour
    % d�terminer les coefficients associ�s.
    du=1;
    try
        [Kc,y0c,ac,t0c]=regression_logistique(t,yc,[],[],[],[],[],[],1,0);
        ycc=fonction_logistique([Kc,y0c,ac],t,t0c);
        if max(abs(ycc-yc))>1
            du=0;
        end
    catch
        du=0;
    end
    if du==0
        Kc=[];
    end
else
    resnormc=[];
end
if isempty(anfin)
    yfin=[];
else
    yfin=feval('fonction_logistique',x,anfin,t0,cho);
end
if grap
    if modmal
        t0b=t(1);
        [A,B]=regression_lineaire(t-t0b,log(y));
        r=A;
        logy0=B;
        y0b=exp(logy0);
        ym=y0b*exp((r*(t-t0b)));
    end
    if isempty(anfin)
        tc=linspace(min(t),max(t),P);
    else
        tc=linspace(min(t),anfin,P);
    end
    hold on;
    verhulst='mod�le de Verhulst';
    if ~isempty(yc)
        if isempty(Kc)
            tb=t;
            yb=yc;
        else
            tb=tc;
            yb=feval('fonction_logistique',[Kc,y0c,ac],tb,t0c,cho);
        end
        if modmal
            plot(t,y,'-*',tb,yb,tc,feval('fonction_logistique',x,tc,t0,cho),t,ym,'-.');
            legend('donn�es','premier mod�le de Verhulst',verhulst,'mod�le de Malthus','location','NorthWest');
        else
            plot(t,y,'-*',tb,yb,tc,feval('fonction_logistique',x,tc,t0,cho));
            legend('donn�es','premier mod�le de Verhulst',verhulst,'location','NorthWest');
        end
    else
        if modmal
            plot(t,y,'-*',tc,feval('fonction_logistique',x,tc,t0),t,ym,'-.');
            legend('donn�es',verhulst,'mod�le de Malthus','location','NorthWest');
        else
            plot(t,y,'-*',tc,feval('fonction_logistique',x,tc,t0));
            legend('donn�es','mod�le','location','NorthWest');
        end
    end
    xlabel('ann�es');
    if isempty(ylabelp)
        ylabel('population');
    else
        ylabel(['population (',ylabelp,')']);
    end
    if tansym
        tc=t0+(log(K/y0-1))/a;
        auxi=a*K/4;
        coeff=[auxi,-1,K/2-auxi*tc];
        trace_droite_fenetre(coeff,'k');
        plot(tc,K/2,'ok');
    end
    hold off;
end