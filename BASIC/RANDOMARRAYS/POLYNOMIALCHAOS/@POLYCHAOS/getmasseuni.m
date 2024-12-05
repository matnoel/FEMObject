function masseuni=getmasseuni(PC,i)

if nargin==1
    
    masseuni=PC.masseuni;
    
elseif length(i)==1
    
    masseuni=PC.masseuni{i};
    
else
    
    masseuni=PC.masseuni(i);
    
end