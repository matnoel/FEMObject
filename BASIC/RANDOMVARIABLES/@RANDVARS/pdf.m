function P = pdf(u,x)

P = zeros(size(x));
for i=1:length(u.RV)
   if isa(u.RV{i},'CONDRANDVAR')
       error('pas prevu pour les conditionnelles')
   else
   P(:,i)=pdf(u.RV{i},x(:,i));
   end
end
P = prod(P,2);


