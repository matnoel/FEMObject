function Z = integrate(a)
% function Z = integrate(a)

LI = a.D*integrate(a.L);
if isa(a.V,'cell')
Z = a.V{1}*LI(1);
for i=2:length(a.V)
Z = Z+  a.V{i}*LI(i);   
end
else
  Z=  double(a.V)*LI;
end

