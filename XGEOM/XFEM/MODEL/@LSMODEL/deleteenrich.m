function S = deleteenrich(S)

S = setfieldgroupelem(S,'lsenrich',0);
g = getnumgroupelemwithfield(S,'lstype','touchcut');
S = setfieldgroupelem(S,'lstype','indomain',g);
g = getnumgroupelemwithfield(S,'lstype','touchbicut');
S = setfieldgroupelem(S,'lstype','indomain',g);
S = concatgroupelem(S);

