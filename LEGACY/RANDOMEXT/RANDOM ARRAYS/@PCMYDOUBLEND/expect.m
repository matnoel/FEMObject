function Z = expect(A,varargin)
% function Z = expect(A,b)
% calcul de Z = E(A*b)
%
% function Z = expect(A,b,c)
% calcul de Z = E(A*b*c)


if ~isa(A,'PCMYDOUBLEND')
    error('pas programme')
end

if isradial(A)
    Ap=PCRADIALMATRIX(A);
else
    Ap=PCMATRIX(A);
end
Ap=reshape(Ap,prod(size(Ap)),1);

Z = expect(Ap,varargin{:});

if isa(Z,'MULTIMATRIX')
    sm = sizem(Z);
    Z = reshape(full(double(Z)),[size(A),prod(sm)]);
    Z = MULTIMYDOUBLEND(Z,length(size(A))+1,sm);
    Z = permutemultidim(Z,A.stodim);
else
    Z = reshape(full(double(Z)),size(A));
    Z = MYDOUBLEND(Z);
end

