function A = calc_ximasse(A,varargin)

A.L = actualise_masse(A.L,varargin{:}); %c'est le getL
A.POLYCHAOS = getPC(A.L);
A.L = calc_ximasse(A.L,varargin{:});

if normest(speye(size(A.D,1))-A.D)>0
    A.DLmasse = multimtimes(A.D,getximasse(A.L)); 
else
    A.DLmasse = getximasse(A.L);
end
