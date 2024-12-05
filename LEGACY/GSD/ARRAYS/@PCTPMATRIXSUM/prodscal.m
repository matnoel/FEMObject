function z = prodscal(x,a)

if isa(x,'PCTPMATRIXSUM') && isa(a,'PCTPMATRIX')
z=0;
for i=1:length(x.funs)
    z = z + prodscal(x.funs{i},a);
end
elseif isa(a,'PCTPMATRIXSUM') && isa(x,'PCTPMATRIX')
z=0;
for i=1:length(a.funs)
    z = z + prodscal(x,a.funs{i});
end
elseif isa(a,'PCTPMATRIXSUM') && isa(x,'PCTPMATRIXSUM')
z=0;
for i=1:length(x.funs)
for j=1:length(a.funs)
    z = z + prodscal(x.funs{i},a.funs{j});
end
end

else
    error('pas programme')
end
    
