function [T,Xc,Yc,er,Tp,erTp,equ]=lotka_volterra_complet(a,b,c,d,x0,y0,tdeb,tfin,seuil,K,te,xe,ye,casrk4,h,xmin,xmax,nx,ymin,ymax,ny,chvit,cylog,traceeq,tracetot,Tp0)

% calcul et tracé des courbes temporelles et dans le plan pour les trois modèle de Lotka_Voltera
%
% *********************************************************
% [T,Xc,Yc,er,Tp,erTp,equ]=lotka_volterra_complet(a,b,c,d,x0,y0,tdeb,tfin,seuil,K,te,xe,ye,casrk4,h,xmin,xmax,nx,ymin,ymax,ny,chvit,cylog,traceeq,tracetot,Tp0) :
%
% Calcul des solutions du modèle complet
% Pour les proies :
%        x'(t) =0                  si x<=seuil
%              =x*a*(1-x/K)        si x>seuil et y<=seuil
%              =x*(a*(1-x/K)-b*y)  si x>seuil et y>seuil
% Pour les prédateurs :
%        y'(t) =0                  si y<=seuil
%              =-y*c               si y>seuil et x<=seuil
%              =+y*(-c+d*x)        si y>seuil et x>seuil
%	* Entrées :
%     * a,b,c,d : paramètres strictement positifs du modèle de
%       Lotka-Voltera de base
%     * x0,y0 : conditions initiales
%	* Entrées optionnelles :
%     * tdeb,tfin : temps de debut et de fin de calcul.
%           * si te est non vide, tdeb=te(1) et tfin=te(end)
%           * si te est vide et seuil=-inf et K=inf (voir
%             ci-dessous) et (x0~=c/d ou y0~=a/b), la période Tp est calculée
%             (sauf si celle-ci est déjà entrée, cf ci-dessous)
%             et 
%                 si tdeb est vide, tdeb=0 et tfin=Tp
%                 sinon tfin=tdeb+Tp
%     * seuil : réel strictement positif ou égal à -inf par défaut (pas de
%       seuil)
%     * K : réel strictement positif ou égal à inf par défaut (pas de capacité
%       de charge)
%     * te,xe,ye : données expérimentales, vides par défaut
%     * casrk4=1 si calcul avec et rk4 et casrk4=0 avec ode45, égal à 1 par
%       défaut. Curieux mais ode45 est moins bon que rk4 !!
%     * pas de discrétisation, égal à 1e-3 par défaut.
%       ATTENTION, peut être légèrement modifié: voir label alpha
%     * xmin,xmax,nx : valeur de x minimal et maximal et nombre de point de
%       discrétisation pour le champ de vitesse.
%       si xmin et xmax non présents, xmin=0
%                                     xmax=max(x);
%                                     nx=40;
%     * ymin,ymax,ny : idem pour y
%     * chvit : entier égal à 0 par défaut,(pas de champ de vitesse) ou 1
%     (champ de vitesse)
%     * cylog : entier égal à 0, par défaut (axes des cycles normaux) ou 1
%     (axes des cycles logarithmique) et si pas de champ de vitesse tracé
%     * traceeq : entier égal à 0, par défaut (pas de tracé de l'équilibre) ou 1
%     (tracé de l'équilibre)
%     * tracetot : entier égal à 1, par défaut (tracé total) ou 0 (pas de tracé)
%     * Tp0 : période de modèle de Lotka-Voltera de base, déjà calculée par
%       la fonction periode_lotka_volterra, vide par défaut
%  * Sorties optionnelles
%     * T,Xc,Yc : temps, effectifs des proies et des prédateurs
%     * er : réel strictement positif : erreur commise par le schéma
%       numérique
%     * Tp, erTp : éventuelles période et erreur (voir fonction
%       periode_lotka_volterra), déterminé dans le cas du modèle de
%       Lotka-Volterra de base (seuil=-inf et K=inf)
%     * equ : les coordonnées de l'équilibre
%
%
% ************ Fonctions auxiliaires utilisées ************
%
%    rk4n
%    periode_lotka_volterra
%
% *********************************************************
%
% (c) 2025 by      Jérôme BASTIEN,
%                  LIBM, Polytech, Université Lyon 1
%                  E-Mail : jerome.bastien@univ-lyon1.fr
%
% Exemples
% a)
% clear all;close all;lotka_volterra_complet(1,1,1,1,0.5,0.5,[],10);
% b)
% clear all;close all;lotka_volterra_complet(1,1,1,1,0.5,0.5);
% c)
% clear all;close all;
% identification_lobry_2017_tot;
% close all;
% [T,Xc,Yc,er]=lotka_volterra_complet(a,b,c,d,x0,y0,[],[],[],[],te,xe,ye);
% d)
% clear all;close all;lotka_volterra_complet(1,1,1,1,0.5,0.5,[],100,[],100);

if nargin<=8||isempty(seuil)
    seuil=-inf;
end
if nargin<=9||isempty(K)
    K=inf;
end
if nargin<=13||isempty(casrk4)
    casrk4=1;
end
if nargin<=14||isempty(h)
    h=1e-3;
end
if nargin<=15
    xmin=[];
end
if nargin<=21||isempty(chvit)
    chvit=0;
end
if nargin<=22||isempty(cylog)
    cylog=0;
end
if nargin<=23||isempty(traceeq)
    traceeq=0;
end
if nargin<=24||isempty(tracetot)
    tracetot=1;
end
if nargin<=25
    Tp0=[];
end
if cylog==1&&chvit==1
    error('On ne peut tracer un champ de vitesse en échelle logarithmique');
end
lvb=~isfinite(seuil)&&~isfinite(K);
% Attention : mettre eps ?!!!
testpart=(x0==c/d&&y0==a/b);
if ~lvb||testpart
    if nargout>=5
        Tp=NaN;
        if nargout>=6
            erTp=NaN;
        end
    end
end
chou=0;
if lvb
    if isempty(Tp0)
        if ~testpart
            [Tp,erTp]=periode_lotka_volterra(a,b,c,d,x0,y0);
        end
    else
        Tp=Tp0;
        erTp=[];
    end
    % label alpha
    if ~testpart
        N=ceil(Tp/h)-1;
    end
end
if nargin>=11&&~isempty(te)
    tfin=te(end);
    tdeb=te(1);
else
    if nargin<=6||isempty(tdeb)
        tdeb=0;
    end
    if nargin<=7||isempty(tfin)
        if lvb&&~testpart
            tfin=tdeb+Tp;
            chou=1;
        else
            error('tfin non défini');
        end
    end
end
if ~lvb
    % label alpha
    N=ceil((tfin-tdeb)/h);
    T=linspace(tdeb,tfin,N);
end
if chvit
    Fx=@(x,y) (x>seuil).*x.*(a*(1-x/K)-b*(y>seuil).*y);
    Fy=@(x,y) (y>seuil).*y.*(-c+d*(x>seuil).*x);
end
F=@(t,Y)[(Y(1)>seuil).*Y(1)*(a*(1-Y(1)/K)-b*(Y(2)>seuil).*Y(2));(Y(2)>seuil).*Y(2)*(-c+d*(Y(1)>seuil).*Y(1))];
if nargout>=7||(traceeq&&tracetot)
    equ=([a/K b;d 0]\[a;c]).';
end
Y0=[x0;y0];
if lvb&&testpart
    if nargin>=11&&~isempty(te)
        Zb=[x0,y0];
        Zb=Zb(1:length(te),:);
    end
    Xc=x0*[1;1];
    Yc=y0*[1,1];
    T=[tdeb;tfin];
    if nargout>=4
        er=0;
    end
else
    chu=0;
    if lvb&&((tfin-tdeb)>=Tp||chou)
        chu=1;
        if chou==1
            qTp=1;
        else
            qTp=floor((tfin-tdeb)/Tp);
        end
        T=linspace(tdeb,tdeb+Tp-h,N);
    end
    if casrk4
        [T,Z]=rk4n(F,[T(1),T(end)],Y0,N);
    else
        %[T,Z]=ode113(F,T,Y0);
        [T,Z]=ode45(F,T,Y0);
        if nargin>=11&&~isempty(te)
            [Tb,Zb]=ode45(F,te,Y0);
        end
    end
    % label alpha
    h0=(T(end)-T(1))/(N-1);
    Xct=Z(:,1);
    Yct=Z(:,2);
    if chu
        Zx=Z(:,1);
        Zx=Zx(:,ones(1,qTp));
        Zx=[Zx(:);Zx(1)];
        Zy=Z(:,2);
        Zy=Zy(:,ones(1,qTp));
        Zy=[Zy(:);Zy(1)];
        Z=[Zx,Zy];
        T=(linspace(tdeb,tdeb+(qTp)*Tp,size(Zx,1))).';
        if T(end)<tfin
            ie=floor((tfin-T(end))/h);
            T=[T;T(end)+(h*(1:ie)).'];
            Z=[Z;Z(1:ie,:)];
        end
        if T(end)<tfin
            if casrk4
                [Tft,Zft]=rk4n(F,[T(end),tfin],Z(end,:),1);
            else
                [Tft,Zft]=ode45(F,[T(end),(T(end)+tfin)/2,tfin],Z(end,:));
                Zft=Zft(end,:);
            end
            T=[T;tfin];
            Z=[Z;Zft];
        end
    end
    if casrk4&&nargin>=11&&~isempty(te)
        Zb=interp1(T,Z,te);
    end
    Xc=Z(:,1);
    Yc=Z(:,2);
    if nargout>=4
        er=erreur(Xct,Yct);
    end
end
% Graphiques
if tracetot
    if chvit&&nargin<=15||isempty(xmin)
        nx=40;
        ny=40;
        xmin=0;
        xmax=max(Xc);
        ymin=0;
        ymax=max(Yc);
    end
    if chvit
        x=linspace(xmin,xmax,nx);
        y=linspace(ymin,ymax,ny);
        [X,Y]=meshgrid(x,y);
    end
    clf;
    if chvit
        hold on;
        quiver(X,Y,Fx(X,Y),Fy(X,Y),'k');
    end
    if cylog==1
        if traceeq
            loglog(Y0(1),Y0(2),'*r',equ(1),equ(2),'or',Xc,Yc);
        else
            loglog(Y0(1),Y0(2),'*r',Xc,Yc);
        end
    else
        if traceeq
            plot(Y0(1),Y0(2),'*r',equ(1),equ(2),'or',Xc,Yc);
        else
            plot(Y0(1),Y0(2),'*r',Xc,Yc);
        end
    end
    if chvit
        if traceeq
            legend('Champ de vitesse','Condition initiale','Equilibre','Trajectoire','location','BestOutside');
        else
            legend('Champ de vitesse','Condition initiale','Trajectoire','location','BestOutside');
        end
    else
        if traceeq
            legend('Condition initiale','Equilibre','Trajectoire','location','BestOutside');
        else
            legend('Condition initiale','Trajectoire','location','BestOutside');
        end
    end
    if cylog==1
        xlabel('effectifs des proies (logarithmiques)');
        ylabel('effectifs des prédateurs (logarithmiques)');
    else
        xlabel('effectifs des proies');
        ylabel('effectifs des prédateurs');
    end
    if chvit
        hold off;
    end
    figure;
    hold on;
    xlim([T(1),T(end)]);
    if nargin>=11&&~isempty(te)
        if traceeq
            plot(T,Xc,'r',T,Yc,'b',te,xe,'r*',te,ye,'bd',[T(1),T(end)],[equ(1),equ(1)],[T(1),T(end)],[equ(2),equ(2)]);
        else
            plot(T,Xc,'r',T,Yc,'b',te,xe,'r*',te,ye,'bd');
        end
    else
        if traceeq
            plot(T,Xc,T,Yc,[T(1),T(end)],[equ(1),equ(1)],[T(1),T(end)],[equ(2),equ(2)]);
        else
            plot(T,Xc,T,Yc);
        end
    end
    if nargin>=11&&~isempty(te)
        teb=[te.';te'];
        xei=[xe.';(Zb(:,1)).'];
        yei=[ye.';(Zb(:,2)).'];
        plot(teb,xei,'r',teb,yei,'b');
    end
    if nargin>=11&&~isempty(te)
        if traceeq
            legend('proies','prédateurs','proies (mesurées)','prédateurs  (mesurés)','proies (équilibre)','prédateurs (équilibre)','location','Best');
        else
            legend('proies','prédateurs','proies (mesurées)','prédateurs  (mesurés)','location','Best');
        end
    else
        if traceeq
            legend('proies','prédateurs','proies (équilibre)','prédateurs (équilibre)','location','Best');
        else
            legend('proies','prédateurs','location','Best');
        end
    end
    xlabel('temps');
    ylabel('effectifs');
    hold off;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fonctions nichées
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function er=erreur(Xc,Yc)

        if casrk4==1
            Xdc=(-Xc(5:end-1)+8*Xc(4:end-2)-8*Xc(2:end-4)+Xc(1:end-5))/(12*h0);
            Ydc=(-Yc(5:end-1)+8*Yc(4:end-2)-8*Yc(2:end-4)+Yc(1:end-5))/(12*h0);
        else
            Xdc=(-Xc(5:end)+8*Xc(4:end-1)-8*Xc(2:end-3)+Xc(1:end-4))/(12*h0);
            Ydc=(-Yc(5:end)+8*Yc(4:end-1)-8*Yc(2:end-3)+Yc(1:end-4))/(12*h0);
        end
        if casrk4==1
            xu=Xc(3:end-3);
            yu=Yc(3:end-3);
        else
            xu=Xc(3:end-2);
            yu=Yc(3:end-2);
        end
        if ~isfinite(seuil)
            er=max([max(abs(xu.*(a*(1-xu/K)-b*yu)-Xdc)),max(abs(yu.*(-c+d*xu)-Ydc))]);
        else
            ind=xu<=seuil;
            if any(ind)
                era=max(abs(Xdc(ind)));
            else
                era=0;
            end
            ind1=xu>seuil&yu<=seuil;
            if any(ind1)
                erb=max(abs(xu(ind1).*(a*(1-xu(ind1)/K))-Xdc(ind1)));
            else
                erb=0;
            end
            ind2=xu>seuil&yu>seuil;
            if any(ind2)
                erc=max(abs(xu(ind2).*(a*(1-xu(ind2)/K)-b*yu(ind2))-Xdc(ind2)));
            else
                erc=0;
            end
            erx=max([era erb erc]);
            ind=yu<=seuil;
            if any(ind)
                era=max(abs(Ydc(ind)));
            else
                era=0;
            end
            if any(ind1)
                erb=max(abs(yu(ind1)*c+Ydc(ind1)));
            else
                erb=0;
            end
            if any(ind2)
                erc=max(abs(yu(ind2).*(-c+d*xu(ind2))-Ydc(ind2)));
            else
                erc=0;
            end
            ery=max([era erb erc]);
            er=max([erx ery]);
        end
    end
end