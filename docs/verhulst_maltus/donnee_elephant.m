% donn�es issues de 
% https://www.mathemathieu.fr/component/attachments/download/294
% https://www.mathemathieu.fr/
% Johan MATHIEU, professeur de math�matiques en coll�ge et lyc�e.
% La dynamique des populations, mod�le continu, Partie 1 : mod�les sansinteraction. TS � D.M.
% ann�es
t=[1905 1923 1930 1939 1945 1950:10:2000];
% population (nombre d'�l�phants)
y=[10 13 29 450 980 3010 5800 6500 7400 7200 7310];
% Attention, ces derniers affich�s � l'unit� pr�s.
% Attention, yc et Kc recalcul�s en interne par regression_logistique, mais
% presque identiques.
% les donn�es ci-dessous conformes au pdf
yc=[10 146 402 1346 2623 3994 6271 7186 7428 7484 7496];
Kc=7500;