function u = getddlname(u,num)
if nargin==1
    num = 1:u.nbddl ;
end
if nargin==2 && numel(num)==1
    u=u.ddl{num};
else
    u=u.ddl(num);
end