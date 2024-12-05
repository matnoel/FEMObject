

m = numel(l);

L0 = expect(A,l,l);
if m>1
    L0 = assembleblock(L0);
end

f0 = expect(bsub,l);
if m>1
    f0 = assembleblock(f0);
end
V = L0\f0;
V = reshape(V,[n,m]);

V = orth(full(V));

fU = V'*bsub;
aU = V'*A*V;

[l,flag]=localstosolverupdate(aU,fU);
l = expand(l);

Rayg = double(full(expectmtimes(l,fU')));
rayg = trace(Rayg);


testsub = 0;
if testsub==1
Vold = V;
lold = l ; 
[RWV,RWD]=eig(Rayg);
V = Vold*RWV;
fU = V'*bsub;
aU = V'*A*V;
[l,flag]=localstosolver(aU,fU);
l = expand(l);
Rayg = real(double(full(expectmtimes(l,fU'))));
Rayg
for i=1:size(V,2)
Ui{i}=V(:,i);
end
for i=1:size(V,2)
fUi{i}=Ui{i}'*bsub;
UAUij{i,i}=Ui{i}'*A*Ui{i};
[li{i},flag]=localstosolver(UAUij{i,i},fUi{i});
CU{i} = expectmtimes(A,li{i},li{i})*Ui{i};
DU{i} = expectmtimes(bsub,li{i});
end

for i=1:size(V,2)
for j=1:size(V,2)
UAU{j,i}=Ui{j}'*A*Ui{i};
UCU(j,i)=Ui{j}'*CU{i};
UDU(j,i)=Ui{j}'*DU{i};
RUU(j,i)=full(expectmtimes(UAU{j,i},li{i},li{j}));
end
RUi(i) = full(expectmtimes(li{i}',fUi{i}'));
end

ZZ = zeros(size(V,2),size(V,2));
for i=1:size(V,2)
for j=1:size(V,2)
ZZ(i,j)=expectmtimes(UAU{i,j},fUi{i},fUi{j});
end
end

end


teststabilityUi=0;
if teststabilityUi
U = V(:,i);
clear RR
for k=1:4
fU=U'*bsub;
aU=U'*A*U;
[l,flag]=localstosolver(aU,fU);
U=expectmtimes(A,l,l)\expectmtimes(bsub,l);
RR(k)=full(expectmtimes(l,fU'))
U = U/norm(U);
end

end

%disp('erreur sur l''erreur relative en utilisant somme(R(Ui))')
%(sum(sum(RUU))-sum(diag(RUU)))/norm(PCRADIALMATRIX(Ui,[n,1],li)-uref,A)^2

