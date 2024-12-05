function rank=getrank(H,node,c,rank)
% Je sais plus trop a quoi ca sert...

% VARIABLES GLOBALES
if nargin==1
    c    = getconnect(H.tree); 
    rank = zeros(1,length(c));
    node = 1;
    rank(node)=H.m;
end

for d=1:H.dim
    subnodes=find(c==node);
    for m=1:H.m
        if isa(H.F{m,d},'HSEP')
            rank = max(rank, getrank( H.F{m,d},subnodes(d),c,rank ) );
        elseif isa(H.F{m,d},'SEP')
            rank(subnodes(d)) = max(rank(subnodes(d)),H.F{m,d}.m);
        end
    end
end