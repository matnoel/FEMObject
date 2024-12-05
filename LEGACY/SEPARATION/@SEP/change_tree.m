function H=change_tree(H,V)
% H=change_tree(H,V)
%     Transforme H du mm type que V
TYPE=class(H);
if isa(H,'HSEP')&&isa(V,'HSEP')
    TYPE=str2func(TYPE);
    H=TYPE(H,V.tree);
elseif isa(H,'HSEP')
    TYPE=TYPE(2:end);
    TYPE=str2func(TYPE);
    H=TYPE(H);
elseif isa(V,'HSEP')
    TYPE=['H' TYPE];
    TYPE=str2func(TYPE);
    H=TYPE(H,V.tree);
else
    % H et V sont des SEP...
end
    