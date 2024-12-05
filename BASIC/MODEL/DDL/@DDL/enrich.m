function u = enrich(u,n,num)
% function u = enrich(u,n,num)


if nargin <=2
    num = 1:u.nbddl ;
end

if nargin==1
    n=1;
end

for i=1:n  
    if n==1
    name = 'E';
    else
    name = ['E' num2str(i)];    
    end

rep = length(u.ddl)+[1:length(num)];
for k=1:length(num)
u.ddl = [u.ddl , {[name u.ddl{num(k)}]}]; 
end
    
u.ddlgroup = [u.ddlgroup , u.ddlgroup(1)];
u.repddlgroup = [u.repddlgroup , {rep}];
u.enrich = [u.enrich,ones(1,length(num))];
end

u.nbddl = length(u.ddl);

