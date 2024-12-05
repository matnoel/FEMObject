function PC = POLYCHAOS(h,p,varargin)
% function PC = POLYCHAOS(h,p)
% h RANDPOLYS ou RANDPOLY ou RANDVARS
% function PC = POLYCHAOS(h,p,'typebase',typebase)
%  typebase=1 : p ordre maxi des polynomes multidimensionnel (par defaut)
%  typebase=2 : p ordre maxi des polynomes unidimensionnel
%       p peut etre un vecteur : ordre different par dimension

if nargin==0;
    PC.indices = zeros(0,1);
    PC.M = 0;
    PC.p = zeros(0,1);
    PC.n = zeros(0,1);
    PC.P = -1;
    PC.typebase = 1;
    PC.masseuni = {};
    PC.masse = {};
    PC.metricuni = {};
    PC.metric = {};
    h=RANDPOLYS();
    PC = class(PC,'POLYCHAOS',h);
    superiorto('RANDPOLYS')
    
else
    
    if nargin==1 && isa(h,'POLYCHAOS')
        PC=h;
    else
        
        if isa(h,'double')
            for kk=1:h
                hh{kk}=POLYHERMITE();
            end
            h = hh;
        end
        
        if ~isa(h,'RANDPOLYS')
            h=RANDPOLYS(h);
        end
        M = getM(h);
        
        typebase=getcharin('typebase',varargin,1);
        % germ=getcharin('germ',varargin);
        
        if nargin==1 || isempty(p)
            p=1;
        end
        if length(p)==1
            p=repmat(p,1,max(1,M));
        end
        
        n = p+1;
        for k=1:M
            if isa(h{k},'POLYFE')
                n(k) = n(k)*getparam(h{k},'n');
            elseif isa(h{k},'POLYLAGRANGE')
                p(k) = getnbpoints(h{k})-1;
                n(k) = p(k)+1;
            elseif isa(h{k},'POLYFELAGRANGE')
                p(k) = getparam(h{k},'m')-1;
                n(k) = getparam(h{k},'n')*getparam(h{k},'m')-(getparam(h{k},'n')-1);
            elseif isa(h{k},'POLYWAVELETS')
                if getparam(h{k},'p')~=p(k)
                    h{k} = setnumber(POLYWAVELETS(getparam(h{k},'n'),p(k)),getnumber(h{k}));
                end
                n(k) = 2^getparam(h{k},'n')*n(k);
            end
        end
        
        if ~ischarin('noindices',varargin)
            PC.indices = PCbase_indices(M,max(n)-1,typebase);
            if length(n)>=M
                for k=1:M
                    rep=find(PC.indices(:,k)>n(k)-1);
                    PC.indices(rep,:)=[];
                end
            else
                error('ordres mal definis')
            end
        else
            PC.indices = zeros(0,M+1);
        end
        
        PC.M = M;
        PC.p = p;
        PC.n = n;
        if ~ischarin('noindices',varargin)
            PC.P = size(PC.indices,1)-1;
        elseif typebase==2
            PC.P = prod(n)-1;
        elseif length(unique(n))==1
            temp = unique(n)-1;
            PC.P = factorial(temp+PC.M)/factorial(temp)/factorial(PC.M)-1;
        else
            PC.P=[];
            warning('pas de calcul de P')
        end
        
        
        PC.typebase = typebase;
        PC.masseuni = {};
        PC.masse = {};
        PC.metricuni = {};
        PC.metric = {};
        % PC.germ = germ;
        
        PC = class(PC,'POLYCHAOS',h);
        superiorto('RANDPOLYS')
        
    end
    
end
