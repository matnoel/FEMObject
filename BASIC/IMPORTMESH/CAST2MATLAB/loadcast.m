function [S,f,cl,mat] = loadcast(file,option)
% function [S,f,cl,mat] = loadcast(file,option)
% -----------      CHARGEMENT de donnees issues de CASTEM    ------------
%
% function [S,f,cl,mat] = loadcast(file)
%    ---> appel de castem si option = 1 : appel de la fonction lance_castem
%                  (castem doit creer les fichiers textes decrits plus loin)
%    ---> appel de cast2matlab_model    : lecture du fichier modele (file_model.txt)
%                                           -->  creation de S.model, S.node, S.elem
%    ---> appel de cast2matlab_load     : lecture des fichiers de chargements  (file_loadi.txt  (i=1...n))
%                                           -->  creation du tableau de cellules f{1} ... f{n}
%    ---> appel de cast2matlab_dis      : lecture des fichiers de conditions cinematiques   (file_disi.txt (i=1...n))
%                                           -->  creation du tableau de cellules cl{1} ... cl{n}
%    ---> appel de calc_S_caract pour remettre en forme et completer le modele
%    ---> appel de calc_cl : interpretation des conditions cinematiques
%
% ---- EXEMPLE COMMENTE DE CREATION DE MODELE SOUS CASTEM :
% ( modele elastique isotrope 2D en contraintes planes )
%  1) donner les options du modele et de la geometrie
%       * OPTI DIME 2 MODE PLAN DEFO TRI3;
%  2)  creer un maillage MAI puis un modele modele MOD
%       * MOD = MODE MAI 'MECANIQUE' 'ELASTIQUE';
%  3)  creer un champ par element DAT contenant toutes les donnees modele et materiau
%       * DAT = MANU CHML MAI 'PLAN' 0 'DEFO' 0 'MATT' 1. 'ELAS' 0. 'ISOT' 0. 'YOUN' 1. 'NU' 0.3 ;
%    Remarques : - il n'est pas necessaire de mettre DIME sachant PLAN
%                - mettre une valeur quelconque pour les champs n'en
%                  necessitant pas (PLAN, TRID, DEFO, CONT, ELAS, ISOT)
%                - MATT indique le numero j du materiau. Ses proprietes
%                  seront dans la cellule mat{j}
%                  On peut associer un materiau par partie de maillage
%                - si la ligne est trop longue, creer DAT en plusieurs etapes
%                - Le type d'element n'est pas toujours necessaire car il
%                  peut etre contenu dans le maillage (sauf pour BARR, POUT, ... qu'il faut alors inserer)
%  4) Sortir le fichier
%       * OPTI SORT 'file_model.txt' ; SORT 'AVS' MAI DAT ;
%  5) Definir des chargements Fi (i=1...n) (champ par point) (a l'aide de FORCE, FSURF, PRES, ...)
%     puis les ecrire :
%       * OPTI SORT 'file_loadi.txt' ; SORT 'AVS' MAI Fi;
%  5) Definir manuellement les conditions cinematiques avec des champs par points DISi.
%     Par exemple :
%       * DIS1 = MANU CHPO 1 L1 'UX' 0 ; bloque le ddl UX sur la partie de maillage L1;
%       * DIS2 = MANU CHPO 2 L2 'UX' 0 'UY' 0 ; bloque le ddl UX sur la partie de maillage L2;
%       * DIS3 = MANU CHPO 3 L3 'DEPL' 0 ; bloque tous les ddl en deplacement sur la partie de maillage L3;
%     Puis sortir les fichiers
%       * OPTI SORT 'file_disi.txt' ; SORT 'AVS' DISi;
%     Remarque : surtout ne pas sortir de maillage ou de modele avec les DISi.

if nargin<2 || option==1
    disp('--------------------------- Appel de Castem ---------------------------')
    lance_castem(file);
    disp('------------------------  Fin appel de Castem -------------------------')
end


S = cast2matlab_model([file '_model.txt']);

if isa(S,'double')
    error(['Le fichier ' file '_model.txt' 'n''existe pas']);
end


f=cell(0,1);
for l=1:inf
    try
        f{l} = cast2matlab_load([file '_load' num2str(l) '.txt']);
    catch
        break
    end
end

cl=cell(0,1);
for l=1:inf
    try
        cl{l} = cast2matlab_dis([file '_dis' num2str(l) '.txt'],S);
    catch
        break
    end
end


%[S,mat] = calc_S_caract(S) ;
%[S,cl] = calc_cl(S,cl);

