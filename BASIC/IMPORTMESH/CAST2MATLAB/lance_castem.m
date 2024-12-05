function lance_castem(file,option)
% function lance_castem(file,option)
% dos(['c:\Calcul\Castem\bin\castemNT .\EXAMPLES\' file],option)
% file peut ne pas avoir l'extension .dgibi



file = strcat('E:\PROGRAMMES\CASTEM\EXAMPLES\',file);

try
    commande = ['c:\Calcul\Castem\bin\castem05 ' file];
    if nargin==1
        option = '-echo';
    end
    dos(commande,option);
catch
    texte=['Le fichier "  ' file '  " n''existe pas'];
    disp(texte); 
end
