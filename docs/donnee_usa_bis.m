% données issues de 
% https://www.imo.universite-paris-saclay.fr/~daniel.perrin/interdisciplines/Cours5equadiff3-2011.pdf
% https://www.imo.universite-paris-saclay.fr/fr/perso/daniel-perrin/
% Daniel Perrin
% Cours numéro 5 :
% Equations différentielles du premier ordre  3, l'équation logistique
%
% Laboratoire de mathématiques d'Orsay
% 
% Université Paris-Saclay
% 91405 Orsay 
% années
t=[1790 1800 1810 1820 1830 1840 1850 1860 1870 1880 1890 1900 1910 ...
1911 1912 1913 1914 1915 1916 1917 1918 ...
1920 1930 1940 1950 1960 1970 1980 1990 2000];
% population en milliers
y=[3929 5308 7240 9638 12866 17069 23192 31443 38558 50156 62948 75995 91972 ...
93512 95055 96599 98144 99690 101235 102779 104320 ...
106022 123203 132165 151326 179323 203302 226342 248710 281422];
y=y/1000; % Pour avoir des millions