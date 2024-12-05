function x = PCTPMATRIXSUM(a,varargin)

if nargin==0
    x.funs = cell(1,0);    
    PC = POLYCHAOSTP();
    x = class(x,'PCTPMATRIXSUM',PC);
    superiorto('PCTPMATRIX')

elseif nargin==1 && isa(a,'POLYCHAOSTP') && ~isa(a,'PCTPMATRIX')
    x.funs = cell(1,0);       
    x = class(x,'PCTPMATRIXSUM',getPC(a));
    superiorto('PCTPMATRIX')
elseif nargin==1 && isa(a,'PCTPMATRIXSUM')
    x=a;

elseif isa(a,'PCTPMATRIX')
    x.funs = cell(1,nargin);
    x.funs{1} = a; 
    PC = getPC(a);
    for i=1:length(varargin)
        if isa(varargin{i},'PCTPMATRIX')  
            x.funs{i+1}=varargin{i};
        else
            error('rentrer des PCTPMATRIX')
        end
    end
    x = class(x,'PCTPMATRIXSUM',PC);
    superiorto('PCTPMATRIX')
else
    error('bad argument')
end
