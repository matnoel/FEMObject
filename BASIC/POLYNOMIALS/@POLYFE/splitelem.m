function h = splitelem(h,e)

param = getparam(h);
I = param.I;
I1 = zeros(0,2);

for i=1:size(I,1)
if ~ismember(i,e) 
    I1 = [I1;I(i,:)];
else
   Ie = I(i,:);
   Ie = linspace(Ie(1),Ie(2),3);
   for k=1:length(Ie)-1
   I1 = [I1;Ie(k),Ie(k+1)];
   end   
end
end
param.I = I1;
param.dx = param.I(:,2)-param.I(:,1);
param.n = length(param.dx);

h = setparam(h,param);