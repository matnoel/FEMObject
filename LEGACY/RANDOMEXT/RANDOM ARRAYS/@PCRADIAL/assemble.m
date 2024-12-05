function K = assemble(Ksto)

fprintf('Assemblage %s',class(Ksto))

mom = getmasse(Ksto) ;
P=getP(Ksto);
n = size(Ksto.V);
m=length(Ksto);


K= spalloc(n(1)*(P+1),n(2)*(P+1),10*max(n)*(P+1)); 
for i=0:P
for j=i:P
pourcentage(i*(P+1)+j+1,(P+1)^2);    
for k=1:m
if abs(mom{k}(i+1,j+1))>eps
K(n(1)*i+1:n(1)*(i+1),n(2)*j+1:n(2)*(j+1))=...
K(n(1)*i+1:n(1)*(i+1),n(2)*j+1:n(2)*(j+1))+...
mom{k}(i+1,j+1)*double(Ksto{k});
end
end
K(n(1)*j+1:n(1)*(j+1),n(2)*i+1:n(2)*(i+1))=...
    K(n(1)*i+1:n(1)*(i+1),n(2)*j+1:n(2)*(j+1));
end
end

