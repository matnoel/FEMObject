function av=moment2(apc)

if iscell(apc)
av = multisum(power(apc.MULTIMATRIX,2));    
else
av = reshape(sum(double(apc).^2,2),size(apc));
end
