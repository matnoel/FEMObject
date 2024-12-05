function n = norm(u)
% function n = norm(u)

n = sqrt(sum(sum(times(u,u),1),2));
