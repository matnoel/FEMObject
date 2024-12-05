function ut = truncate(u,modes)
% function ut = truncate(u,modes)
% format of modes:   modes={1:10,4:8,3:6};
%       If modes={1:10,[],3:6}, then there is
%       no truncature in the dim 2

dim=u.dim;

emptymodes=cellfun(@isempty,modes);
for i=1:dim
    if emptymodes(i)
        modes{i}=1:u.m(i);
    end
end

F=cell(dim,1);
for i=1:dim
    F{i}=u.F{i}(modes{i});
end

switch dim
    case 2
        alpha=u.alpha(modes{1},modes{2});
    case 3
        alpha=u.alpha(modes{1},modes{2},modes{3});
    otherwise
        error('not implemented')
end



ut=TSEPMATRIX(F,tensor(alpha,cellfun(@numel,modes)));
