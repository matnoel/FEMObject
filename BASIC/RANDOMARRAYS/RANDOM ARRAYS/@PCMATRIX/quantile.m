function a = quantile(apc,q,n)
% function a = quantile(apc,q,n)
% apc : PCMATRIX
% q : quantiles
% n : nombre de tirages pour le calcul des quantiles


a = random(apc,n);

if isa(a,'MULTIMATRIX') & iscell(a)
a = cell2mat(a); 
elseif ~isa(a,'double') && ~isa(a,'MULTIMATRIX')
error('mauvais type')
end

a = quantile(double(a)',q)';
    
    
if length(q)>1
a = MULTIMATRIX(a,size(apc),[length(q),1]);
else
a = reshape(a,size(apc));    
end