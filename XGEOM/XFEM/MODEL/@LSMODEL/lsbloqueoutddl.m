function S=lsbloqueoutddl(S)

elemout = getnumgroupelemwithfield(S,'lstype','out');
elemnonout = setdiff(1:getnbgroupelem(S),elemout);

numelemnonout = getnumelem(S,elemnonout,'local');

connec = getconnec(S);
nodenonout = find(sum(connec.node2elem(numelemnonout,:),1));
nodeout = setdiff(1:getnbnode(S),nodenonout);

[ddlbloque,ddltemp] = findddl(S,'all',nodeout);

S.MODEL = addbc(S.MODEL,ddlbloque,zeros(length(ddlbloque),1),ddltemp);
