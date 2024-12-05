function av = variance(apc)
% function av = variance(apc)

if iscell(apc)
    av = multisum(power(apc.MULTIMATRIX,2));
    av = av-mean(apc).^2;
else
    av = reshape(sum(double(apc).^2,2),size(apc));
    av = av-mean(apc).^2;
end
