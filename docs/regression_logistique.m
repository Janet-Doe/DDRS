function [K,y0,a,t0,resnorm,resnormc,yfin]=regression_logistique(t,y,modmal,P,yc,ylabelp,MaxFunEvals,cho,iter,grap,anfin,tansym)

% paramètres de K,y0,a,t0 de la regression logistique (modèle de Verhulst) pour des données (t,y)
%
% *********************************************************
%   [K,y0,a,t0,resnorm,resnormc,yfin]=regression_logistique(t,y,modmal,P,yc,ylabelp,MaxFunEvals,cho,iter,grap,anfin,tansym)
%
%	* Entrées :
%     * t,y : deux tableaux ligne représentant les données ti,yi (dates,
%       population)
%	* Entrées optionnelles :
%     * modmal entier dans {0,1} égal à 1 par défaut
%      * si 0 : pas d'affichage du modèle de Malthus
%      * si 1 : affichage du modèle de Malthus
%     * P : nombre de points pour le graphique continu, égal à 1000 par
%     défaut
%     * yc : autres mesures faites (avec les mêmes ordonnées t, qui donnera
%       l'erreur en sortie resnormc), vide par défaut.
%     * labelp : coefficient multiplicateur pour la légende "population" en
%     y, vide par défaut.
%     * MaxFunEvals: vide par défaut, pour l'option de lsqcurvefit
%     * cho : voir fonction fonction_logistique
%     * iter : entier dans {0,1}, égal à 0 par défaut.
%        affiche les résultats partiels de l'optimisation si iter=1.
%     * grap : entier dans {0,1}, égal à 1 par défaut.
%        produit le graphique si grap=1. Sur ce graphique, apparaissent :
%          * les données mesurées (ti,yi)
%          * la courbe correpondant au modèle de Verhulst (voir la fonction
%            fonction_logistique)
%          * éventuellement, la courbe correpondant au modèle de Malthus
%          * la courbe correpondant au modèle de Verhulst  dont les valeurs
%             sont dans yc
%     * anfin : année de fin de graphique, qui n'est pas nécessairement la
%       donnée maximale de t. vide par défaut. Si non vide, la valeur de la
%       population pour cette année est renvoyée dans yfin
%     * tansym : entier dans {0,1}, égal à 0 par défaut.
%        si =1, affiche la tangente au point d'inflexion.
%  * Sorties
%     * K,y0,a,t0 tel que yi=K/(1+(K/y0-1)*exp(-a*(t-t0)) au sens des
%       moindres carrés, correspondant à la fonction logistique (voir
%       fonction_logistique)
%     Attention, t0 est arbitraire et sera choisi égal à la première valeur de t
%  * Sorties optionnelles
%     * resnorm : sommes des carrés des écarts entre les données et le
%       modèle.
%     * resnormc : voir yc
%     * yfin : voir anfin
%
%
% ************ Fonctions auxiliaires utilisées ************
%
%    fonction_logistique
%    regression_lineaire
%    trace_droite_fenetre
%
% *********************************************************
%
% (c) 2024 by      Jérôme BASTIEN,
%                  LIBM, Polytech, Université Lyon 1
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
    error('Il existe des valeurs négatives ou nulles de y');
end
% Initialisation avec la Méthode de Perrin
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
    % déterminer les coefficients associés.
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
    verhulst='modèle de Verhulst';
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
            legend('données','premier modèle de Verhulst',verhulst,'modèle de Malthus','location','NorthWest');
        else
            plot(t,y,'-*',tb,yb,tc,feval('fonction_logistique',x,tc,t0,cho));
            legend('données','premier modèle de Verhulst',verhulst,'location','NorthWest');
        end
    else
        if modmal
            plot(t,y,'-*',tc,feval('fonction_logistique',x,tc,t0),t,ym,'-.');
            legend('données',verhulst,'modèle de Malthus','location','NorthWest');
        else
            plot(t,y,'-*',tc,feval('fonction_logistique',x,tc,t0));
            legend('données','modèle','location','NorthWest');
        end
    end
    xlabel('années');
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