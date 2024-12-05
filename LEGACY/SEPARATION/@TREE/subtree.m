function sub_T = subtree(T,i)
% sub_T = sous-arbre de l'arbre T a partir du noeud i
%
%             T
%           / | \
%          /  |  \
%         /   |   \i
%        /|  /\   sub_T

%%% Trouver les noeuds en dessous
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

    if test == 1
        pass = [pass ii];
    end
end

%%% Reconstruire les liens
sub_connect = zeros(1,size(pass,2));
for ii=1:size(pass,2)
    if pass(ii) == i
        link=0;
    else
        oldlink = T.connect(pass(ii));
        link = find(pass==oldlink);
    end
    sub_connect(ii)=link;%% Son nouveau lien
end



sub_T = TREE(sub_connect);




