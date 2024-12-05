function M = create_groupddl(M)

[node,nbddl] = create_groupddlenrich(getnode(M),calcgroupddl(M),getlevelsets(M));

M = setnode(M,node);
M = setnbddl(M,nbddl);
