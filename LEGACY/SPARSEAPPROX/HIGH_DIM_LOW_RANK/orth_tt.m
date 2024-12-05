function v = orth_tt(v)
% function v = orth_tt(v)

d = numel(v);

% last core
sz = size(v{d});
w = reshape(v{d},sz(1),sz(3));
[q,r] = qr(w',0);
sz(1) = size(r,1);
v{d} = reshape(q',sz);
v{d-1} = ttm(v{d-1},r,2);

% other cores
for i = d-1:-1:2
    sz = size(v{i});
    vv = reshape(v{i},sz(1),sz(2)*sz(3));
    [q,r] = qr(vv',0);
    sz(1) = size(r,1);
    v{i} = reshape(q',sz);
    v{i-1} = ttm(v{i-1},r,2);
end

end
