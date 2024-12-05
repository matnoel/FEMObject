function PCTP = POLYCHAOSTP(h,p,varargin)
%function PCTP = POLYCHAOSTP(h,p,'groups',groups)
% create a model for tensor product of stochastic function spaces
% groups: cell array containing the groups of stochastic dimensions
% groups{i} is the ith set of stochastic dimensions
% examples: 
% - POLYCHAOSTP(h,p,'groups',{[1,2],[3,4,5],[6:n]}) creates a tensor product
% of three spaces V_1\otimes V_2\otimes V_3, with V_1 associated with
% dimensions [1,2], V_2 with dimensions [3,4,5] ... approximation spaces
% are associated with multidimensional hermite polynomials with partial
% degree p
% - POLYCHAOSTP(h,p,'groups',{[1,2],[3,4,5],[6:n]},'typebase',1) is the same except
% that each space V_i are associated with multidimensional hermite
% polynomials with total degree p
% - POLYCHAOSTP(RANDPOLYS(POLYLEGENDRE(),5),p,{[1,2],[3,4,5]},'typebase',1) 
% ...


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

    if length(p)==1
        p=repmat(p,1,max(1,m));
    end

    PC = POLYCHAOS(h,p,'typebase',typebase,'noindices');
    h = RANDPOLYS(PC);
    PCTP.PCgroups = cell(1,length(PCTP.groups));
    for k=1:length(PCTP.groups)
        g = PCTP.groups{k};
        PCTP.PCgroups{k} = POLYCHAOS(h(g),p(g),'typebase',typebase);    
    end
    PCTP = class(PCTP,'POLYCHAOSTP',PC);
    superiorto('POLYCHAOS')
    PCTP = calc_mean(PCTP);
end

