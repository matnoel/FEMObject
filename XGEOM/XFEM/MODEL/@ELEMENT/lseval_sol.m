function [qP,repP] = lseval_sol(elem,node,q,P,ddl,ls)

[repP,numelem]=ispointin(elem,node,P);

P=P(repP);

if isempty(repP)
    qP = zeros(length(DDL(ddl)),size(q,2),0);
else
    
    
    elem = getelem(elem,numelem,'global','nounique');
    xnode = getcoord(node,getconnec(elem)');
    
    qe=localize(elem,q);
    
    xi = calc_xi(elem,xnode,P);
    
    if isenrich(elem)
        ls=getlevelset(ls,getlsnumber(elem));
        N = calc_Nls(elem,xnode,xi,ls);
    else
        N = calc_N(elem,xnode,xi);
    end
    
    
    qP = double(N*qe);
    
    if nargin>=5
        ddlnode = getddlnode(elem);
        repddl=findddl(ddlnode,ddl);
        qP = qP(repddl,:,:);
    end
    
end