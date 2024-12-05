function x = getnbddlfree(S)
% function x = getnbddlfree(S)

if ~isempty(S.bc)
   switch S.bc.type 
       case 'dirichlet'
           x = getnbnode(S)-length(S.bc.node);
       case 'periodic'
           x = getnbnode(S)-1;
   end
else
    x = getnbnode(S);
end