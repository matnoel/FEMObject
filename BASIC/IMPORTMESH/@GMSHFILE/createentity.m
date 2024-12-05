function u = createentity(u,name,values,number)
% function u = createentity(u,name,values,number)

if nargin==4
    u = stepcounter(u,number);
else
    u = stepcounter(u);
    number = u.counter;
end

s = entity(name,number,values);
u = addstring(u,s);
