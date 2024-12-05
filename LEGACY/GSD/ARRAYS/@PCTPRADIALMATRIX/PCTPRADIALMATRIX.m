function R = PCTPRADIALMATRIX(varargin)


if nargin==0
    R.V = MULTIMATRIX();
    R.L = PCTPMATRIX();
    R.m=[];
    R.DLmasse = cell(1,0);
    R.D=[];    

    R= class(R,'PCTPRADIALMATRIX',POLYCHAOSTP());
    superiorto('PCTPMATRIX','PCTPMATRIXSUM')
elseif nargin==1 && isa(varargin{1},'PCTPRADIALMATRIX')
    R = varargin{1};
elseif nargin==1 && isa(varargin{1},'PCTPMATRIX')
    R = varargin{1};
    R = PCTPRADIALMATRIX(getphi0(R),size(R),setphi0(R,1));
elseif nargin==1 && isa(varargin{1},'PCTPMATRIXSUM')
    R = varargin{1};
    funs = getfuns(R);
    L = cell(1,length(funs));
    V = cell(1,length(funs));
    for i=1:length(funs)
        L{i} = setphi0(funs{i},1);
        V{i} = getphi0(funs{i});
    end
    R = PCTPRADIALMATRIX(V,size(R),L);
elseif nargin==2
    sV=varargin{1};
    PC = POLYCHAOSTP(varargin{2});
    R.V = MULTIMATRIX(sparse(prod(sV),0),sV,[0,1]);
    R.L = cell(0,1);
    R.m=0;
    R.DLmasse = cell(1,0);
    R.D=speye(0,0);    

    R= class(R,'PCTPRADIALMATRIX',PC);
    superiorto('PCTPMATRIX','PCTPMATRIXSUM')
elseif nargin>=2
    V=varargin{1};
    if isa(V,'cell') || isa(V,'double')    
        sV=varargin{2};
        V=MULTIMATRIX(V,sV);
    elseif ~isa(V,'MULTIMATRIX')
        error('mauvais argument')
    end
    L = varargin{3};

    if  ~isa(L,'cell') 
        if isa(L,'PCTPMATRIX') || isa(L,'PCTPMATRIXSUM')
            LL = cell(1,numel(L));
            for i=1:numel(LL)
                LL{i} = L(i);
            end
            L=LL;
        else
            error('L must be a PCTPMATRIX OU PCTPMATRIXSUM or a cell containing them')
        end
    end
    L=L(:);

    if length(V)~=size(L,1)
        error('number of L and V functions must equal')      
    end

    R.V = V;
    R.L = L;
    R.m = length(V);
    R.V = reshapem(R.V,[R.m,1]);
    PC = getPC(R.L{1});
    R.DLmasse=cell(1,0);
    if nargin==4
        R.D = varargin{4};
    else
        R.D=speye(R.m);
    end

    R= class(R,'PCTPRADIALMATRIX',PC);
    superiorto('PCTPMATRIX','PCTPMATRIXSUM')

else
    error('bad arguments')

end

