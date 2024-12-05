function A = calc_masse(A,varargin)

A.L = actualise_masse(A.L,varargin{:});
A.POLYCHAOS = getPC(A.L);

A.L = calc_ximasse(A.L,varargin{:});
A.DLmasse = multimtimes(A.D,getximasse(A.L)); 
