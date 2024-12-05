function v=addddl(v,u,varargin)

if isa(u,'DDLVECT') || isa(u,'DDLSCAL') || isa(u,'DDLTENS2')  || isa(u,'DDLTENS4')
    v.ddl=[v.ddl,getddl(u)];
    v.ddlgroup = [v.ddlgroup,{u}];
    v.repddlgroup =  [v.repddlgroup,{v.nbddl+(1:getnbddl(u))}];
elseif isa(u,'char')
    v.ddl=[v.ddl,{u}];
    v.ddlgroup = [v.ddlgroup,{'char'}];
    v.repddlgroup =  [v.repddlgroup,{v.nbddl+1}];
elseif isa(u,'DDL')
    v.ddl=[v.ddl,u.ddl];
    v.ddlgroup = [v.ddlgroup,{u}];
    v.repddlgroup =  [v.repddlgroup,{v.nbddl+(1:getnbddl(u))}];
else
    error('pas prevu')
end

v.nbddl = length(v.ddl);

if nargin>=3
    v = addddl(v,varargin{:});
end
