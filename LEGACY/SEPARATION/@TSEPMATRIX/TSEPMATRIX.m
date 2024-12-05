function A = TSEPMATRIX(F,alpha,varargin)
% function A = TSEPMATRIX(F,alpha)
% function A = TSEPMATRIX(dim,m)


if nargin==0
    A.dim = 0;
    A.m = 0;
    A.alpha = tensor([]);
    A.F = {};
    A = class(A,'TSEPMATRIX');
    
    
elseif nargin==1 && isa(F,'TSEPMATRIX')
    A = F;
    
    
elseif nargin==1 && isa(F,'struct')
    A = class(F,'TSEPMATRIX');
    
    
elseif isa(F,'double')
    A.dim = F;
    if nargin==1
        alpha=0;
        A.m = zeros(1,A.dim);
    else
        A.m=alpha;
    end
    if A.dim==1
        A.alpha = tensor(zeros(1,alpha),alpha);
    else
        A.alpha = tensor(zeros(alpha),alpha);
    end
    A.F = cell(A.dim,1);
    for i=1:A.dim
        A.F{i}=cell(A.m(i),1);
    end
    A = class(A,'TSEPMATRIX');
    
    
elseif isa(F,'cell')
    if nargin~=2
        dim=size(F,1);
        m=zeros(1,dim);
        for i=1:dim
            m(i)=size(F{i},1);
        end
        alpha=tenones(m);
    end
    A=TSEPMATRIX(ndims(alpha));
    A.F=F;
    A.alpha=alpha;
    A.m=size(alpha);
    
elseif isa(F,'ttensor')
    A=TSEPMATRIX();
    A.alpha=F.core;
    A.F=F.U';
    A.m=size(A.alpha);
    A.dim=ndims(A.alpha);
    A=splitvectors(A);
    
elseif isa(F,'SEPMATRIX')
    A.dim=F.dim;
    A.m=F.m*ones(1,F.dim);
    A.alpha=sptensordiag(F.alpha,F.dim);
    A.F=cell(A.dim,1);
    for i=1:A.dim
        A.F{i}=cell(F.m,1);
        for j=1:F.m
            A.F{i}{j}=F.F{j,i};
        end
    end
    A = class(A,'TSEPMATRIX');
    
elseif isa(F,'ttensor')
    A.dim=ndims(F);
    A.m=F.core.size;
    A.alpha=F.core;
    
    A.F=cell(A.dim,1);
    for i=1:A.dim
        A.F{i}=cell(A.m(i),1);
        for j=1:A.m(i)
            A.F{i}{j}=F.U{i}(:,j);
        end
    end
    A = class(A,'TSEPMATRIX');
    
end


superiorto('SEPSOLVER')

