function metricuni=getmetricuni(PC,i)

if nargin==1
    
    metricuni=PC.metricuni;
    
elseif length(i)==1
    
    metricuni=PC.metricuni{i};
    
else
    
    metricuni=PC.metricuni(i);
    
end