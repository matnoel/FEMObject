function rad=truncate(rad,k)
%function rad=truncate(rad,k)
if length(k)==1
    scank=1:k;
else
    scank=k(:)';
end
k=length(scank);
DL = rad.D * rad.L ; 
L = DL(scank);

if isa(rad.V,'cell')
V=cell(1,k);

%L=cell(1,k);
for i=1:k
   V{i}=rad.V{scank(i)}; 
%   L{i}=DL(scank(i));
end
else
V = double(rad.V);
V = V(:,scank);
V = MULTIMATRIX(V,size(rad.V));
end

rad = TIMERADIALMATRIX(V,size(rad),L);
   