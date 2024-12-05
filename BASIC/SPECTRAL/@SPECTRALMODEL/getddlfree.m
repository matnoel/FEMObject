function x = getddlfree(S)
% function x = getnbddlfree(S)
if ~isempty(S.bc)
   switch S.bc.type 
       case 'dirichlet'
           x = setdiff(1:getnbnode(S),S.bc.node);
       case 'periodic'
           x = 1:getnbnode(S)-1;
   end
else
    x = 1:getnbnode(S);
end