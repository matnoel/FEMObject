function C = times3D(A,B)
% function C = times3D(A,B)

sA = [size2D(A),sizeND(A)];
sB = [size2D(B),sizeND(B)];

if sA(3)==1
    C = A*B;
elseif sA(2)==1 && sA(1)==1
    C = prod3Ddoublemat(A,B);
elseif sB(2)==1 && sB(1)==1
    C = prod3Ddoublemat(B,A);
elseif sA(2)==sB(1)
    C = prod3Dmatmat(A,B);
else
    error('la multiplication n''est pas possible')
end

function C = prod3Ddoublemat(A,B)
sA = size3D(A);
sB = size3D(B);
B = reshape(B,sB(1)*sB(2),sB(3));
A = repmat(reshape(A,1,sB(3)),[sB(1)*sB(2),1]);
C = reshape(A.*B,sB);

return

function C = prod3Dmatmat(A,B)
sA = size3D(A);
sB = size3D(B);

% A = repmat(A,[sB(2),1,1]);
% A = permute(A,[1,3,2]);
% A = reshape(A,sA(1)*sB(2)*sA(3),sA(2));
%
% B = permute(B,[2,3,1]);
% B = reshape(B,[1,sB(2),sB(1)*sB(3)]);
% B = repmat(B,[sA(1),1,1]);
% B = reshape(B,sA(1)*sB(2)*sB(3),sB(1));
% C = sum(A.*B,2) ;
% C = reshape(C,[sA(1,sB(2),sA(3)]);

C = zeros(sA(1),sB(2),sA(3));
A = permute(A,[2 1 3]);
for i=1:sA(1)
    for j=1:sB(2)
        C(i,j,:) = sum(A(:,i,:).*B(:,j,:),1);
    end
end
return
