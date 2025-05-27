function trace_multi_cycle(Xc,Yc,equ,cylog,ch)

% Trace de plusieurs cycles :
%
% *********************************************************
% trace_multi_cycle(Xc,Yc,equ,cylog,nomx,nomy) :
%
%	* Entrées :
%     * Xc,YC : cellules contenant les abscisses et ordonnées (chaque
%     cycles est donné par Xc{i},Yc{i}, pas nécessairement de la même
%     taille)
%	* Entrées optionnelles :
%     * equ : les coordonnées de l'équilibre, vide ou inexistant par défaut
%       si de taille (1,2) : le même équilibre pour tout
%       sinon de taille (Nt,2) où Nt est le nombre de courbe.
%     * cylog : entier égal à 0, par défaut (axes des cycles normaux) ou 1
%     * ch : chaîne pour les différentes légendes, vide par défaut
%
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
if nargin<=2
    equ=[];
end
if nargin<=3||isempty(cylog)
    cylog=0;
end
if nargin<=4
    ch=[];
end
Nt=length(Xc);
hold on;
if ~isempty(equ)
    q=size(equ,1);
else
    q=0;
end
colors=hsv(Nt);
if q==0
    for i=1:Nt
        if cylog==1
            loglog(Xc{i},Yc{i},'color',colors(i,:));
        else
            plot(Xc{i},Yc{i},'color',colors(i,:));
        end
    end
    for i=1:Nt
        auxix=Xc{i};
        auxiy=Yc{i};
        if cylog==1
            loglog(auxix(1),auxiy(1),'color',colors(i,:),'Marker','*');
        else
            plot(auxix(1),auxiy(1),'color',colors(i,:),'Marker','*');
        end
    end
    if ~isempty(ch)
        chi=cell(1,Nt);
        for i=1:Nt
            chi{i}=ch{i};
        end
        legend(chi,'location','BestOutside');
    end
elseif q==1
    if cylog==1
        loglog(equ(1),equ(2),'ok');
    else
        plot(equ(1),equ(2),'ok');
    end
    for i=1:Nt
        if cylog==1
            loglog(Xc{i},Yc{i},'color',colors(i,:));
        else
            plot(Xc{i},Yc{i},'color',colors(i,:));
        end
    end
    for i=1:Nt
        auxix=Xc{i};
        auxiy=Yc{i};
        if cylog==1
            loglog(auxix(1),auxiy(1),'color',colors(i,:),'Marker','*');
        else
            plot(auxix(1),auxiy(1),'color',colors(i,:),'Marker','*');
        end
    end
    if ~isempty(ch)
        chi=cell(1,Nt+1);
        chi{1}='Equilibre commun';
        for i=1:Nt
            chi{i+1}=ch{i};
        end
        legend(chi,'location','BestOutside');
    end
else
    for i=1:Nt
        if cylog==1
            loglog(Xc{i},Yc{i},'color',colors(i,:));
        else
            plot(Xc{i},Yc{i},'color',colors(i,:));
        end
    end
    for i=1:q
        if cylog==1
            loglog(equ(i,1),equ(i,2),'o','color',colors(i,:));
        else
            plot(equ(i,1),equ(i,2),'o','color',colors(i,:));
        end
    end
    for i=1:Nt
        auxix=Xc{i};
        auxiy=Yc{i};
        if cylog==1
            loglog(auxix(1),auxiy(1),'color',colors(i,:),'Marker','*');
        else
            plot(auxix(1),auxiy(1),'color',colors(i,:),'Marker','*');
        end
    end
    if ~isempty(ch)
        chi=cell(1,Nt);
        for i=1:Nt
            chi{i}=ch{i};
        end
        legend(chi,'location','BestOutside');
    end
end
xlabel('effectifs des proies');
ylabel('effectifs des prédateurs');
hold off;