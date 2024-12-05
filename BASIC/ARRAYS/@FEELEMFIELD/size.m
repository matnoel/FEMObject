function siz = size(u,n)
% function siz = size(u,n)

for k=1:length(u.value)
    siz{k} = [size(u.value{k},1),size(u.value{k},2)];
    if nargin==2
        siz{k} = siz{k}(n);
    end
end
ok = 1;sizall = siz{1};
for k=2:length(u.value)
    ok = ok & all(siz{k}==sizall);
end
if ok
    siz = sizall;
else
    error('les size associees aux differents groupes d''elements ne correspondent pas')
end
