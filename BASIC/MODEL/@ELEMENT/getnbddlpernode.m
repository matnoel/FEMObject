function n = getnbddlpernode(elem,type)
% function n = getnbddlpernode(elem,type)

if nargin==1
    type = 'ddlnode';
end
ddlchoice = eval(['@get' type]);

n = length(ddlchoice(elem));
