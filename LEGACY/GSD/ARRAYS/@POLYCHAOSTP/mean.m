function Hm = mean(PC,i)

if nargin==1
Hm = PC.mean ; 
elseif length(i)==1  
Hm = PC.mean{i};
else
Hm = PC.mean(i);
end
