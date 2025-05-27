function lotka_volterra_bis(a,b,c,d,x0,y0,tdeb,tfin,seuil,K,te,xe,ye,eqcomp,chmax)

% calcul et trac� des courbes temporelles et dans le plan pour les trois mod�les de Lotka_Voltera pour �ventuellement plusieurs valeurs de param�tres
%
% lotka_volterra_bis(a,b,c,d,x0,y0,tdeb,tfin,seuil,K,te,xe,ye,eqcomp,chmax) :
%
% Calcul des solutions du mod�le complet et affichage des courbes
% temporelles (si une seule valeur de chaque param�tre, voir plus bas, et
% ensemble des cycles).
% Pour les proies :
%        x'(t) =0                  si x<=seuil
%              =x*a*(1-x/K)        si x>seuil et y<=seuil
%              =x*(a*(1-x/K)-b*y)  si x>seuil et y>seuil
% Pour les pr�dateurs :
%        y'(t) =0                  si y<=seuil
%              =-y*c               si y>seuil et x<=seuil
%              =+y*(-c+d*x)        si y>seuil et x>seuil
%	* Entr�es :
%     * a,b,c,d : param�tres strictement positifs du mod�le de
%       Lotka-Voltera de base
%     * x0,y0 : conditions initiales
%     Attention, plusieurs de ces 6 param�tres peuvent �tre un tableaux,
%     tous ces tableaux d�vant �tre de m�me taille.
%	* Entr�es optionnelles :
%     * tdeb,tfin : temps de debut et de fin de calcul.
%           * si te est non vide, tdeb=te(1) et tfin=te(end)
%           * si te est vide et seuil=-inf et K=inf (voir
%             ci-dessous) et (x0~=c/d ou y0~=a/b), la p�riode Tp est calcul�e
%             et
%                 si tdeb est vide, tdeb=0 et tfin=Tp
%                 sinon tfin=tdeb+Tp
%     * seuil : r�el strictement positif ou �gal � -inf par d�faut (pas de
%       seuil)
%     * K : r�el strictement positif ou �gal � inf par d�faut (pas de capacit�
%       de charge)
%      Attention, plusieurs de ces 2 param�tres peuvent �tre un tableaux,
%      tous ces tableaux d�vant �tre de m�me taille.
%     * te,xe,ye : donn�es exp�rimentales, vides par d�faut
%     * eqcomp : cha�ne de caract�re contenant une �quation liant un des
%     param�tre a, b, c, d, x0, y0, seuil, K aux autres
%     * chmax : nombre maximum de courbes pour l'affichage de la l�gende,
%     �gal � 15 par d�faut.
%
%
% ************ Fonctions auxiliaires utilis�es ************
%
% lotka_volterra_complet
% trace_multi_cycle
%
% *********************************************************
%
% (c) 2025 by      J�r�me BASTIEN,
%                  LIBM, Polytech, Universit� Lyon 1
%                  E-Mail : jerome.bastien@univ-lyon1.fr
%
% Exemples
% a)
% clear all;close all;lotka_volterra_bis(1,1,1,1,0.5,0.5,[],10);
% b)
% clear all;close all;lotka_volterra_bis(1,1,1,1,0.5,0.5);
% c)
% clear all;close all;
% identification_lobry_2017_tot;
% close all;
% lotka_volterra_bis(a,b,c,d,x0,y0,[],[],[],[],te,xe,ye);
% d)
% clear all;close all;lotka_volterra_bis(1,1,1,1,0.5,0.5,[],100,[],100);
% e)
% clear all;close all;lotka_volterra_bis(1*linspace(1,2,10),1,1,1,0.5,0.5);
% f)
% clear all;close all;lotka_volterra_bis(1*linspace(1,2,10),1*linspace(1,2,10),1,1,0.5,0.5);
% g)
% clear all;close all;lotka_volterra_bis(1*linspace(1,2,10),1*linspace(1,2,10),1,1,0.5,0.5,[],[],[],[],[],[],[],'d=b');

if nargin<=6
    tdeb=[];
end
if nargin<=7
    tfin=[];
end
if nargin<=8||isempty(seuil)
    seuil=-inf;
end
if nargin<=9||isempty(K)
    K=inf;
end
if nargin<=10
    te=[];
    xe=[];
    ye=[];
end
if nargin<=13
    eqcomp=[];
end
if nargin<=14
    chmax=15;
end
u=8;
N=zeros(1,u);
N(1)=length(a);
N(2)=length(b);
N(3)=length(c);
N(4)=length(d);
N(5)=length(x0);
N(6)=length(y0);
N(7)=length(seuil);
N(8)=length(K);
i=find(N>1);
if isempty(i)
    p=1;
else
    if length(i)>=2
        if max(abs(N(i)-N(i(1))))~=0
            error('plus de deux param�tres varient avec des cardinaux diff�rents');
        end
    end
    p=N(i(1));
end
auxi=ones(1,p);
nom={'a','b','c','d','x0','y0','s','K'};
if N(1)==1
    at=a(auxi);
else
    at=a;
end
if N(2)==1
    bt=b(auxi);
else
    bt=b;
end
if N(3)==1
    ct=c(auxi);
else
    ct=c;
end
if N(4)==1
    dt=d(auxi);
else
    dt=d;
end
if N(5)==1
    x0t=x0(auxi);
else
    x0t=x0;
end
if N(6)==1
    y0t=y0(auxi);
else
    y0t=y0;
end
if N(7)==1
    seuilt=seuil(auxi);
else
    seuilt=seuil;
end
if N(8)==1
    Kt=K(auxi);
else
    Kt=K;
end
Xt=[at;bt;ct;dt;x0t;y0t;seuilt;Kt];
if p==1
    [Tl,Xcl,Ycl,er,Tp,erTp]=lotka_volterra_complet(a,b,c,d,x0,y0,tdeb,tfin,seuil,K,te,xe,ye,[],[],[],[],[],[],[],[],[],[],1,1);
    disp('erreur sch�ma');
    disp(er);
    if ~isnan(erTp)
        disp('�ventuelle erreur p�riode');
        disp(erTp);
    end
else
    Xc=cell(1,p);
    Yc=Xc;
    equ=zeros(p,2);
    for q=1:p
        disp(['Calcul de la courbe num�ro ',int2str(q),' en cours ...']);
        a=at(q);
        b=bt(q);
        c=ct(q);
        d=dt(q);
        x0=x0t(q);
        y0=y0t(q);
        seuil=seuilt(q);
        K=Kt(q);
        if ~isempty(eqcomp)
            eval([eqcomp,';']);
        end
        [Tl,Xcl,Ycl,er,Tp,erTp,equl]=lotka_volterra_complet(a,b,c,d,x0,y0,tdeb,tfin,seuil,K,...
            [],[],[],[],[],[],[],[],[],[],[],[],[],0,0);
        disp('erreur sch�ma');
        disp(er);
        if ~isnan(erTp)
            disp('�ventuelle erreur p�riode');
            disp(erTp);
        end
        Xc{q}=Xcl;
        Yc{q}=Ycl;
        equ(q,:)=equl;
    end
    if max(abs(diff(equ,1)))==0
        equ=equ(1,:);
    end
    if p<=chmax
        ch=cell(1,p);
        for q=1:p
            if i==1
                ch{q}=[nom{i},'=',num2str(Xt(i,q))];
            else
                auxi=[];
                for j=1:length(i)
                    auxi=[auxi,nom{i(j)},'=',num2str(Xt(i(j),q))];
                    if j<=length(i)-1
                        auxi=[auxi,', '];
                    end
                end
                ch{q}=auxi;
            end
        end
    else
        ch=[];
    end
    trace_multi_cycle(Xc,Yc,equ,[],ch);
end