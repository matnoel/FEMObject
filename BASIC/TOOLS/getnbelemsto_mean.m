function [nbelemsto,elemsto] = getnbelemsto_mean(M,n)
% function [nbelemsto]= getnbelemsto_mean(M,n)
% --> nbelemsto : moyenne du nombre de cellules des partitions
% stochastiques d'un groupe d'�l�ments
% M : MODEL al�atoire
% n : n� du groupe --> M.groupelem{n}
% sortie : temp : nombre de cellules pour chaque �l�ment fini du groupe
if israndom(M)
    D = M.groupelem{n};
    elemsto = zeros(getnbelem(D),1);
    for j=1:getnbelem(D)
        stomesh = getparam(getelem(D,j),'stomesh');
        elemsto(j) = stomesh{1}.n;
    end
    nbelemsto = sum(elemsto)/size(elemsto,1);
else
    nbelemsto = 0;
end
