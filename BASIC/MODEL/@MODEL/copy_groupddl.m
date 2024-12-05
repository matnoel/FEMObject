function M = copy_groupddl(M,N)

if isa(N,'NODE')
    [M.node,M.nbddl]=copy_groupddl(M.node,N);
    
elseif isa(N,'MODEL')
    M = copy_groupddl(M,getnode(N));
else
    error('mauvais argument')
end
