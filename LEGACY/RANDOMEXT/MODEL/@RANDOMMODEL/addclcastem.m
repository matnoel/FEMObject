function M = addclcastem(M,cl)
% function M = addclcastem(M,cl)

if isa(cl,'struct')
    cl = {cl};
end

totalbloque = [];
for j=1:M.nbcl
    totalbloque = union(totalbloque,M.cl{j}.ddlbloque);
end

for k=1:length(cl);
    bloque = [];
    ud = zeros(M.nbddl,1);
    
    for l=1:length(cl{k}.ddl)
        switch cl{k}.ddl{l}
            case 'DEPL'
                ddl = DDL(DDLVECT('U',M.syscoord,'TRANS'));
            case 'ROTA'
                ddl = DDL(DDLVECT('R',M.syscoord,'ROTA'));
            otherwise
                ddl = cl{k}.ddl;
        end
        
        for j=1:length(cl{k}.node)
            numddl = findddl(M.groupddl{M.repnodeingroupddl(cl{k}.node(j),1)}.ddlnode,ddl);
            ddlbloque = M.groupddl{M.repnodeingroupddl(cl{k}.node(j),1)}.numddl(M.repnodeingroupddl(cl{k}.node(j),2),numddl) ;
            bloque = union(bloque,ddlbloque);
            ud(ddlbloque) = cl{k}.value(j,:);
        end
    end
    
    ud(totalbloque) = 0; % pour eviter les redondances de conditions cinematiques
    
    M.ddlbloque = union(M.ddlbloque,bloque);
    M.ddlfree = setdiff(1:M.nbddl,M.ddlbloque);
    M.cl{M.nbcl+1}.cl = cl;
    M.cl{M.nbcl+1}.ud = SPACEFUN(ud,'primal');
    M.cl{M.nbcl+1}.ddlbloque = setdiff(bloque,totalbloque);
    M.nbcl = M.nbcl+1;
    
end

% pour avoir des condiditons cinematiques non nulles, il suffit de
% construire une SPACEFUN 'primal' avec les cl{k}.value
