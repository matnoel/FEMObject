function P = inv_r1(P)
% function Q = inv_r1(P)
% invert a rank-1 operator

if P.m>1
    error('P must be a rank-1 operator')
end

P.alpha=1/P.alpha;
for i=1:P.dim
    P.F{i}=inv(P.F{i});
end

