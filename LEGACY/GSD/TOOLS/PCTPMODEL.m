function [X,PC] = PCTPMODEL(a,varargin)
% function [X,PC] = PCTPMODEL(v,'order',p,'groups',groups,'typebase',typebase)
% v : RANDVARS or RANDPOLYS that characterize the underlying stochastic function
% space
% denoting a = RANDVARS(v) (the underlying random variables)
% and h = RANDPOLYS(v) the associated orthonormal polynomials
% (other approximation bases can be provided giving for v a RANDPOLYS)
%
% the following function is called (see POLYCHAOSTP)
% PC = POLYCHAOSTP(RANDPOLYS(a),p,'groups',groups,'typebase',typebase)
% X : contains the representation of the random variables a in the tensor
% product approximation space.
%
%
% different POLYCHAOS can also be provided
% function PC = PCTPMODEL(pc1,pc2,...)
% pci : POLYCHAOS

if isa(a,'RANDVARS') || isa(a,'RANDVAR') || isa(a,'RANDPOLYS')
    p = getcharin('order',varargin);
    if isempty(p)
        error('preciser l''ordre de la decomposition (degre des polynomes)')  
    elseif length(p)==1
        p=repmat(p,1,max(1,getM(a)));
    end
    varargin = setcharin('order',varargin,p);

    if isa(a,'RANDPOLYS')
        H = a;
        a = RANDVARS(H);
    else
        H = create_randpolys(a,varargin{:});
    end

    PC = POLYCHAOSTP(H,p,varargin{:});
    if ~ischarin('nomasse',varargin)
        PC = calc_masseuni(PC,[],'cell2mat');
        PC = calc_massegroups(PC);
        PC = calc_metricgroups(PC);
    end

    X = cell(1,getM(a));
    for i=1:getM(a)
        g = findgroupwithdim(PC,i);
        apc = double(project(a{i},getpcgroup(PC,g)));
        apc = full(apc(:));
        X{i} = full(one(PC));
        X{i} = setphi(X{i},apc,g);
        if ~ischarin('nomasse',varargin) 
            X{i} = calc_ximasse(X{i});
        end
    end

elseif isa(a,'PCMATRIX')


elseif isa(a,'POLYCHAOS')


else


end
