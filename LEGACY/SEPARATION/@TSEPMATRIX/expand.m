function A = expand(u)

gu=gathervectors(u);
switch u.dim
    case 1
        A=gu.F{1}*double(u.alpha);
    case 2
        A=gu.F{1}*double(u.alpha)*gu.F{2}';
    otherwise
        A=ttm(u.alpha,gu.F{1},1);
        for i=2:u.dim
            A=ttm(A,gu.F{i},i);
        end
end
A=double(full(A));

