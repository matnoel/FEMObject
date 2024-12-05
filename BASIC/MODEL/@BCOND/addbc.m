function B = addbc(B,numddl,value,ddl)
% function B = addbc(B,numddl,value,ddl)

B.nbcl = B.nbcl + 1;

[ddlbloque,rep] = setdiff(numddl,B.ddlbloque); % pour eviter les redondances de conditions
B.BC{B.nbcl}.ddlbloque = ddlbloque;
B.BC{B.nbcl}.value = value(rep,:);
B.BC{B.nbcl}.type = 0;
if isa(B.BC{B.nbcl}.value,'double') && all(B.BC{B.nbcl}.value==0)
    B.BC{B.nbcl}.ishomogeneous = 1;
else
    B.BC{B.nbcl}.ishomogeneous = 0;
end
B.BC{B.nbcl}.ddl = ddl;

B.ddlbloque = union(B.ddlbloque,ddlbloque);
B.ddlfree = setdiff(B.ddlfree,B.ddlbloque);
