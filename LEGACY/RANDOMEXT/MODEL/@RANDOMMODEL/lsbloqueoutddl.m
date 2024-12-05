function S=lsbloqueoutddl(S)

nodecut=[];
nodein=[];
nodeout=[];

for p=1:S.nbgroupelem
    elem = S.groupelem{p};
    switch getlstype(elem)
        case {'in','indomain'}
            nodein = union(nodein,unique(getconnec(elem)));
        case 'out'
            nodeout = union(nodeout,unique(getconnec(elem)));
        case 'cut'
            nodecut = union(nodecut,unique(getconnec(elem)));
        otherwise
            error('bad lstype')
    end
    
end

nodeout = setdiff(nodeout,nodecut);
nodeout = setdiff(nodeout,nodein);

ddlbloque = findddl(S,'all',nodeout);
S.BCOND = addbc(S.BCOND,ddlbloque,zeros(length(ddlbloque),1));

