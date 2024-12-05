function M = cast2matlab_load(M,file)
% function M = cast2matlab_load(M,file)

[field,node,elem,nb] = cast2matlab_field(file);

M = addloadcastem(M,field);
