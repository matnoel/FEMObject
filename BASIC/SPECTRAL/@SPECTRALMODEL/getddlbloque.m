function x = getddlbloque(S)
% function x = getnbddlbloque(S)
if ~isempty(S.bc)
   switch S.bc.type 
       case 'dirichlet'
           x = S.bc.node;
       case 'periodic'
           x = getnbnode(S);
   end
else
    x = [];
end