function [N,numnode]=getnodenextto(N,P,choix)
% function [N,numnode]=getnodenextto(N,P,choix)
% choix = 'global' : numero du noeud
% choix = 'local' : position dans la liste de noeuds

if isa(P,'POINT')
    [PN,repnode]=getpointnextto(N.POINT,P);
    numnode = getnumber(N,repnode);  
    N = getnode(N,repnode,'local'); 
    
else
        error('pas programme')
end

if nargin==3 && strcmp(choix,'local')
   numnode = repnode;
end