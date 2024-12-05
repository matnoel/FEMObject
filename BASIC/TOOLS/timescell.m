function C=timescell(A,B)
if ~isa(A,'cell') | ~isa(B,'cell')
    error('A et B doivent etre des cell')
end

sA = size(A) ;
sB = size(B) ;

if prod(sA)~=prod(sB)
    error('cell dimension must agree')    
end


sA1 = size(A{1}) ;
sB1 = size(B{1}) ;

A=reshape([A{:}],[sA1 prod(sA)]);
B=reshape([B{:}],[sB1 prod(sB)]);

C = times3D(A,B);
C = reshape(num2cell(C,[1 2]),sA);


