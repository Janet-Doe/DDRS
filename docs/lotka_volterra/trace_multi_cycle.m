function trace_multi_cycle(Xc,Yc,equ,cylog,ch)

% Trace de plusieurs cycles :
%
% *********************************************************
% trace_multi_cycle(Xc,Yc,equ,cylog,nomx,nomy) :
%
%	* Entr�es :
%     * Xc,YC : cellules contenant les abscisses et ordonn�es (chaque
%     cycles est donn� par Xc{i},Yc{i}, pas n�cessairement de la m�me
%     taille)
%	* Entr�es optionnelles :
%     * equ : les coordonn�es de l'�quilibre, vide ou inexistant par d�faut
%       si de taille (1,2) : le m�me �quilibre pour tout
%       sinon de taille (Nt,2) o� Nt est le nombre de courbe.
%     * cylog : entier �gal � 0, par d�faut (axes des cycles normaux) ou 1
%     * ch : cha�ne pour les diff�rentes l�gendes, vide par d�faut
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
ylabel('effectifs des pr�dateurs');
hold off;