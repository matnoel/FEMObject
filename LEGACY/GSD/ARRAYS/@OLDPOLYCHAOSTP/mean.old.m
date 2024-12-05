function Hm = mean(PC,i)

if nargin==1
Hm = cell(1,getnbdim(PC));
for i=1:getnbdim(PC)
  Hm{i} = mean(getpoly(PC,i),0:getn(PC,i)-1);
end
elseif length(i)==1  
Hm = mean(getpoly(PC,i),0:getn(PC,i)-1);
else
Hm = cell(1,length(i));
for k=1:length(i)
  Hm{k} = mean(getpoly(PC,i(k)),0:getn(PC,i(k))-1);    
end 
end
