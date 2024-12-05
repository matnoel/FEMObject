function PCTP = POLYCHAOSTP(h,p,varargin)

if nargin==0
PCTP = struct();
PCTP.mean = [];
PCTP.groups = {};
PCTP.PCgroups = {};
PCTP = class(PCTP,'POLYCHAOSTP',POLYCHAOS());
superiorto('POLYCHAOS')
elseif nargin==1 && isa(h,'POLYCHAOSTP')
PCTP = h;    
else
PCTP = struct();
PCTP.mean = [];

typebase=getcharin('typebase',varargin,1);
if isa(h,'RANDPOLYS')
 m=getM(h);   
else
 m=h;   
end
PCTP.groups = getcharin('groups',varargin,{mat2cell(1:m,1,ones(1,m))});    
if ~isa(PCTP.groups,'cell')
    error('error dans groups')
end
PC = POLYCHAOS(h,p,'typebase',typebase,'noindices');
h = RANDPOLYS(PC);
PCTP.PCgroups = cell(1,length(PCTP.groups));
for k=1:length(PCTP.groups)
g = PCTP.groups{k};
PCTP.PCgroups{k} = POLYCHAOS(h(g),p(g),'typebase',typebase);    
end
PC = POLYCHAOS(h,p,'typebase',typebase,'noindices');
PCTP = class(PCTP,'POLYCHAOSTP',PC);
superiorto('POLYCHAOS')
PCTP = calc_mean(PCTP);
end

