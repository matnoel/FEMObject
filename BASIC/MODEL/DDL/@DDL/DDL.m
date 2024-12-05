function v = DDL(varargin)
% function v = DDL(DDLSCAL1,DDLVECT1,...)
% arguments en entree : objets de type DDLSCAL ou DDLVECT ou DDLTENS2 ou DDLTENS4

if nargin==0
    v.nbddl = 0;
    v.ddl = cell(1,0);
    v.ddlgroup = cell(1,0);
    v.repddlgroup = cell(0,1);
    v.enrich = zeros(1,v.nbddl);
    v = class(v,'DDL');
elseif nargin==1 && isa(varargin{1},'DDL')
    v = varargin{1};
else
    v = DDL();
    v = addddl(v,varargin{:});
end

% for l=1:nargin
%     if isa(varargin{l},'DDLVECT') || isa(varargin{l},'DDLSCAL') || isa(varargin{l},'DDLTENS2') || isa(varargin{l},'DDLTENS4')
%         v.ddl = [v.ddl,getddl(varargin{l})];
%         v.ddlgroup = [v.ddlgroup,{varargin{l}}];
%         v.repddlgroup =  [v.repddlgroup,{v.nbddl+[1:getnbddl(varargin{l})]}];
%     elseif isa(varargin{l},'char')
%         v.ddl = [v.ddl,{varargin{l}}];
%         v.ddlgroup = [v.ddlgroup,{'char'}];
%         v.repddlgroup =  [v.repddlgroup,{v.nbddl+1}];
%     end
%     v.nbddl = length(v.ddl);
% end
