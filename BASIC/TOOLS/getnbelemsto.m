function [nbelemsto,elemsto] = getnbelemsto(M)
% function [nbelemsto,temp] = getnbelemsto(M)
% --> donne le nombre de cellules stochastiques totale d'un maillage
% �l�ments finis
% M : MODEL al�atoire
% sorties : nbelemsto
%           elemsto : nombre de cellules pour chaque �l�ments finis par
%           groupe d'�l�ments.
if israndom(M)
    nbelemsto = zeros(M.nbgroupelem,1);

    for p=1:M.nbgroupelem
        temp = zeros(getnbelem(M.groupelem{p}),1);
        for j=1:getnbelem(M.groupelem{p})
            stomesh = getparam(getelem(M.groupelem{p},j),'stomesh');
            temp(j) = stomesh{1}.n;
        end
        nbelemsto(p) = sum(temp);
        elemsto{p} = temp; 
    end
    nbelemsto = sum(nbelemsto);
else
    nbelemsto = 0;
end
