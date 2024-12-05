function A=speye(varargin)
PC=POLYCHAOS(getclassin('POLYCHAOS',varargin));
s = getclassin('double',varargin);
A = speye(s);
A = A(:) * mean(PC)';

A = PCMATRIX(A,[s,s],PC);