function CT = cuttree(T,i)
% function CT = cuttree(T,node)
%   -> Coupe l'arbre au noeud 'node' :
leaf=[];
pass=[];
for ii=1:size(T.connect,2)
    sup = (ii);
    test= 0;
    while sup ~= 0
        if sup == i
            test = 1;
        end
        sup= T.connect(sup);
    end
    if test == 0 || ii==i
        pass = [pass ii];
    end
    if test==1 && isempty(find( T.connect==ii , 1))
        leaf = [leaf ii];
    end
end

% Reconstruire les liens
sub_connect = zeros(1,length(pass)+length(leaf));
for ii=1:size(pass,2)
    oldlink = T.connect(pass(ii));
    if oldlink==0
        link = 0;
    else
        link = find(pass==oldlink);
    end
    sub_connect(ii)=link;%% Son nouveau lien
end


% A la fin, on complete avec les feuilles :
NI=find(pass==i);
for ii=1:length(leaf)
    sub_connect(length(pass)+ii)=NI;
end

CT = TREE(sub_connect);



