function y=fonction_logistique(x,xdata,t0,cho)

% fonction logistique y=K/(1+(K/y0-1)*exp(-a*(xdata-t0)) où x=[K,y0,a], vectorielle en xdata
% (cho : argument optionnel, voir dans la fonction)
%
% y=fonction_logistique(x,xdata,t0,cho)

if nargin<=3||isempty(cho)
    cho=1;
end
% test : on vérifie que0<y0<K 
if 0>=x(2)||x(2)>=x(1)
    error('On n''a pas 0<y0<K');
end
% test : on vérifie que 0<a 
if 0>=x(3)
    error('On n''a pas 0<a');
end
y=x(1)./(1+(x(1)/x(2)-1)*exp(-x(3)*(xdata-t0)));
% ce dernier test (y<K) n'est pris en compte que si cho==1
% le test y>0 est inutile car c'est vrai si le test ci-dessus est vrai
if cho==1
    ind=y>=x(1);
    if any(ind)
        error('Il existe des valeurs de y supérieures ou égales à K');
    end
end