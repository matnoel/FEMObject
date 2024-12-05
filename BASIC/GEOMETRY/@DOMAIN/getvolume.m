function v = getvolume(D)
% function v = getvolume(D)

P = getvertices(D);
N = zeros(length(P),1);
for i=1:length(P)
    N(i) = norm(P{i});
end
min_indice = N==min(min(N));
max_indice = N==max(max(N));

bb = P{max_indice}-P{min_indice};
v = 1;
for i=1:length(bb)
    v = v*bb(i);
end
