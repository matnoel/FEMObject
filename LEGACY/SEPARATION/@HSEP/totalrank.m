function TotRank = totalrank(H)
% TotRank est le rang de l'object SEP(H)

connect = getconnect(H.tree);
node    = find(connect==0);
subnode = find(connect==node);
TotRank = 0;
% Parcours des rang (somme)
for r=1:H.m
    subrankr = 1;
    for d=1:size(subnode,2) % Parcours des dimensions (produit)
        if isa(H.F{r,d},'HSEP')
            subrankr = subrankr * totalrank(H.F{r,d});
        else
            subrankr = subrankr * H.F{r,d}.m;
        end
    end
    TotRank = TotRank + subrankr;
end

