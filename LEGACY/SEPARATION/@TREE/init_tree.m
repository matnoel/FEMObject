function [dim,var,var2dim,leaf,Cvar] = init_tree(T)
% Le but est de deduire les attribus dim, var, var2dim, leaf et Cvar a
% partir de connect.
% ATTENTION A LA REDACTION DE CONNECT !!!!!!!!!!!!!!!!!!!!!
% 1 EST TOUJOURS LE PREMIER NOEUD
% UNE BRANCHE SEULE EST TOUJOURS CONSTITUEE DE DEUX NOEUDS...

connect=T.connect;
nbnode = length(connect);

dim  = 0;
leaf = [];

%%%%%%% D'abort la dim et les feuilles
for i=1:nbnode
    subnodes = find(connect == i, 1);
    if isempty(subnodes) %% On est sur une feuille
        % Ajouter une dimension
        dim = dim+1;
        % Retenir que ce neoud est une feuille
        leaf = [leaf i];
    end
end

var     = zeros(nbnode,1);
var2dim = cell(dim,1);
%%%%%%% 
root=find(connect==0);
for i=1:dim
    node   = leaf(i);
    var(node)=i;
    var2dim{i}=[];
    while node ~= root
        subdim = find( find(connect==connect(node))==node);
        var2dim{i} = [subdim var2dim{i}];
        node = connect(node);
    end
end


Cvar  = cell(nbnode,1);
for i=1:length(leaf)
    node   = leaf(i);
    dimnum = var(node);
    while node ~= 0
        Cvar{node} = sort([Cvar{node} dimnum]);
        node = connect(node);
    end
end



