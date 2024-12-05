function u = PCMYDOUBLEND(a,varargin)
% function u = PCMYDOUBLEND(a,PC,stodim)
% a : MYDOUBLEND
% PC : POLYCHAOS
% stodim : indique la dimension stochastique dans a


if nargin==0

    u.stodim = 5;
    u.storage = 'PC';
    u.V = zerosND(0,1,1,1,0);     
    u.L = POLYCHAOS();

    u=class(u,'PCMYDOUBLEND');
    superiorto('PCMATRIX','PCRADIALMATRIX','MYDOUBLEND');

elseif nargin==3 

    u.stodim = varargin{2};

    if length(u.stodim)>1
        error('la dimension stochastique doit etre un entier')
    end
    if isa(a,'MULTIMATRIX')
        a = MYDOUBLEND(a,u.stodim);
    elseif ~isa(a,'MYDOUBLEND')  
        error('rentrer un MYDOUBLEND ou une MULTIMATRIX')
    end    

    s = size(a);
    if u.stodim>length(s)
        s=[s,ones(1,u.stodim-length(s))];
    end
    L = varargin{1};
    if isa(L,'PCMATRIX')
        n = numel(L);
    elseif strcmp(class(L),'POLYCHAOS')
        n = length(L);
    else
        keyboard
    end
%keyboard
    if n~=s(u.stodim)
        error(['le dimension ' num2str(u.stodim) ' du MYDOUBLEND doit correspondre � celle de ' class(varargin{1})])
    end
    if isa(varargin{1},'PCMATRIX')
        u.storage = 'RADIAL';
    elseif strcmp(class(varargin{1}),'POLYCHAOS')
        u.storage = 'PC';
    else
        error('rentrer un POLYCHAOS ou une PCMATRIX comme second argument')
    end
    u.V = a;
    u.L = L ;
    u=class(u,'PCMYDOUBLEND');
    superiorto('PCMATRIX','PCRADIALMATRIX','MYDOUBLEND');
else
    error('mauvais argument')
end

