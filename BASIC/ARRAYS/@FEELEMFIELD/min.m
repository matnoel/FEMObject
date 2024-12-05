function [v,IJ] = min(u)
% function [v,IJ] = min(u)

if ~israndom(u)
    for k=1:length(u.value)
        v{k} = min(u.value{k},[],4);
        [v{k},I{k}] = min(v{k},[],3);
        v{k} = double(v{k});
        siz{k} = size(v{k});
        v{k} = v{k}(:);
        I{k} = I{k}(:);
    end
    ok = 1;sizall = siz{1};
    for k=2:length(v)
        ok = ok & all(siz{k}==sizall);
    end
    if ok
        [v,J] = min([v{:}],[],2);
        v = reshape(v,sizall);
        IJ = zeros(size(J));
        for k=1:length(u.value)
            rep = find(J==k);
            IJ(rep) = I{k}(rep);
        end
        v = reshape(v,sizall);
        IJ = [J,IJ];
    end
else
    v = decompmatrix(getPC(u.value{1}),[],[],@fun,u);
end


function m = fun(a,u)

ur = randomeval(u,a);

m = min(min(ur.value{1},[],3),[],4);
for k=2:length(u.value)
    m = [m, min(min(ur.value{k},[],3),[],4)];
end
m = min(m,[],2);

return
