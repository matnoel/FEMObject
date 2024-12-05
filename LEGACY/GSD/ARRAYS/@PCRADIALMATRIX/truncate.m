function rad=truncate(rad,k)
%function rad=truncate(rad,k)
if length(k)==1
    scank=1:k;
else
    scank=k(:)';
end
k=length(scank);
DL = rad.D * rad.L ; 
V = rad.V{scank};
L = DL(scank);
%for i=1:k
%   V{i}=rad.V{scank(i)}; 
%   L{i}=DL(scank(i));
%end
rad = PCRADIALMATRIX(V,size(rad),L);

   