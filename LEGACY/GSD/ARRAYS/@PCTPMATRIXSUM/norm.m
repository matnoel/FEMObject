function n = norm(x)
% function n = norm(x)

n = 0;

for k1=1:getm(x)
    for k2=1:getm(x)
        n = n + sum(sum(expecttimes(x.funs{k1},x.funs{k2})));
    end
end

n = sqrt(abs(n));
