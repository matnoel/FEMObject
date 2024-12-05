function x = PCMATRIX(a,varargin)
% function x = PCMATRIX(a,s,PC)
% PC : POLYCHAOS
% s taille de la matrice
% a : double de taille prod(s)*length(PC)
if nargin==0
    value = MULTIMATRIX();
    PC = POLYCHAOS();
    x=struct();
    x.ximasse = {};
    x=class(x,'PCMATRIX',value,PC);
elseif isa(a,'PCMATRIX')
    x = a;
elseif isa(a,'MULTIMATRIX')

    if nargin==2 && isa(varargin{1},'POLYCHAOS') 
        PC = getPC(varargin{1});
        s = size(a);
    elseif nargin==3 
        PC = getPC(varargin{2});
        s=varargin{1};
    end

    if length(a)~=length(PC)
        error('la MULTIMATRIX doit etre adapt�e au chaos')
    end    
    if ~isempty(s)
        value = reshape(a,s);
    end
    x=struct();   
    x.ximasse = {};
    x=class(x,'PCMATRIX',value,PC);

elseif nargin==2 && isa(a,'double') && isa(varargin{1},'POLYCHAOS')
    PC = getPC(varargin{1});  
    x = PCMATRIX(sparse(prod(a),length(PC)),a,PC);

else
    x=struct();   
    x.ximasse = {};
    PC=getPC(varargin{2});
    s = varargin{1};
    if ~isa(s,'double') || isempty(s) || numel(s)==1
        error('donner la taille de la matrice')
    end

    if isa(a,'double')
        value = MULTIMATRIX(reshape(a,prod(s),length(PC)),s,[length(PC),1]);
    elseif isa(a,'cell') 
        value = MULTIMATRIX(a,s,[length(PC),1]); 
    else
        error('pas prevu')
    end
    x=class(x,'PCMATRIX',value,PC);
end

