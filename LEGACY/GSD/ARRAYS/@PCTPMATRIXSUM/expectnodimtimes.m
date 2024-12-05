function x = expectnodimtimes(nodim,x,a,b)
% function x = expectnodimtimes(nodim,x,a,b)

if nargin==3
if isa(x,'PCTPMATRIXSUM') && isa(a,'PCTPMATRIX')
for i=1:length(x.funs)
    x.funs{i} = expectnodimtimes(nodim,x.funs{i},a);
end
elseif isa(a,'PCTPMATRIXSUM') && isa(x,'PCTPMATRIX')
for i=1:length(a.funs)
    a.funs{i} = expectnodimtimes(nodim,x,a.funs{i});
end
x=a;

elseif isa(a,'PCTPMATRIXSUM') && isa(x,'PCTPMATRIXSUM')
z = PCTPMATRIXSUM(a.POLYCHAOSTP);
z.funs = cell(1,length(a.funs)*length(x.funs));
ll=0;
for i=1:length(x.funs)
    for j=1:length(a.funs)
        ll=ll+1;
    z.funs{ll} = expectnodimtimes(nodim,x.funs{i},a.funs{j});
    end
end
x=z;

else
    error('pas programme')
end
elseif nargin==4
error('pas programme')    
else
    error('pas programme')
end
