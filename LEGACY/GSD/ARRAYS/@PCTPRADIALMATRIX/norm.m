function n = norm(a,varargin)
% function n = norm(a,varargin)

n = 0;

for k1=1:getm(a)
    for k2=1:getm(a)
        n = n + prodscal(getV(a,k1),getV(a,k2))*prodscal(getL(a,k1),getL(a,k2));
    end
end

n = sqrt(abs(n));
