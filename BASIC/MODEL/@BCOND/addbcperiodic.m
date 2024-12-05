function B = addbcperiodic(B,numddl1,numddl2,ddl)
% function B = addbcperiodic(B,numddl1,numddl2,ddl)

B.nbcl = B.nbcl + 1;

[numddl2,rep] = setdiff(numddl2,B.ddlbloque); % pour eviter les redondances de conditions
numddl1=numddl1(rep);
[numddl1,rep] = setdiff(numddl1,B.ddlbloque); % pour eviter les redondances de conditions
numddl2=numddl2(rep);
B.BC{B.nbcl}.ddlfree = numddl1; 
B.BC{B.nbcl}.ddlbloque = numddl2; 
B.BC{B.nbcl}.type = 1;
B.BC{B.nbcl}.ishomogeneous = 1;
B.BC{B.nbcl}.ddl = ddl;

B.ddlbloque = union(B.ddlbloque,numddl2);
B.ddlfree = setdiff(B.ddlfree,B.ddlbloque);
