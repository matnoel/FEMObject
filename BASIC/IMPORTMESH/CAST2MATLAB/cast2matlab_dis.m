function M = cast2matlab_dis(M,file)
% function M = cast2matlab_dis(M,file)

nbgroup = 0;
[field,node,elem,nb] = cast2matlab_field(file);

cl.ddl = field.chponame;
cl.node = zeros(size(node,1),1);
cl.value = zeros(size(node,1),length(cl.ddl));
[P,rep1,rep2] = intersect(POINT(node),M.node);
cl.node = rep2 ;
cl.value = field.chpo ;

M = addclcastem(M,cl);
