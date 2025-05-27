% données issues de 
% https://www.mathemathieu.fr/component/attachments/download/294
% https://www.mathemathieu.fr/
% Johan MATHIEU, professeur de mathématiques en collège et lycée.
% La dynamique des populations, modèle continu, Partie 1 : modèles sansinteraction. TS – D.M.
% années
t=[1905 1923 1930 1939 1945 1950:10:2000];
% population (nombre d'éléphants)
y=[10 13 29 450 980 3010 5800 6500 7400 7200 7310];
% Attention, ces derniers affichés à l'unité près.
% Attention, yc et Kc recalculés en interne par regression_logistique, mais
% presque identiques.
% les données ci-dessous conformes au pdf
yc=[10 146 402 1346 2623 3994 6271 7186 7428 7484 7496];
Kc=7500;