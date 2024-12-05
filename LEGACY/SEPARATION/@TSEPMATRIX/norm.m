function n = norm(t)
% function n = norm(t)

t = orth(t);
n = norm(t.alpha);
