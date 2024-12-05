function p = RANDPOLY(type,param,domain)

if nargin==0
    p.type = '';
    p.paramnames = {};
    p.param = struct();
    p.domain = [];
    p.number = [];
    p = class(p,'RANDPOLY');    
elseif isa(type,'RANDPOLY')
    p = type ;
else
    p.type = type;
    p.paramnames = fieldnames(param);
    p.param = param;
    p.domain = domain;
    p.number = [];
    p = class(p,'RANDPOLY');
end


