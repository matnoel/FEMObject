function [qP,repP] = eval_sol(elem,node,q,P,ddl)
% function [qP,repP] = eval_sol(elem,node,q,P,ddl)

[repP,numelem]=ispointin(elem,node,P);

P = P(repP);

if isempty(repP)
    qP = zeros(length(DDL(ddl)),size(q,2),0);
else
    elem = getelem(elem,numelem,'global','nounique');
    xnode = getcoord(node,getconnec(elem)');
    
    qe = localize(elem,q);
    
    xi = calc_xi(elem,xnode,P);
    
    N = calc_N(elem,xnode,xi);
    
    qP = double(N*qe);
    
    if islocal(elem)
        R = calc_P(getddlnode(elem),getsyscoord(elem));
        qP = R'*qP;
    end
    
    if nargin>=5
        ddlnode = getddlnode(elem);
        repddl = findddl(ddlnode,ddl);
        qP = qP(repddl,:,:);
    end
    
end