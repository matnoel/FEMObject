function u = getddl(u,num)
% function u = getddl(u,num)

if nargin>1
    u.ddl = u.ddl(num);
    u.nbddl = length(u.ddl);
    u.enrich = u.enrich(num);
end
