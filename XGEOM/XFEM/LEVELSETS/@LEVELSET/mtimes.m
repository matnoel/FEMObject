function ls=mtimes(u,v)
if isa(u,'LEVELSET') && isa(v,'double')
ls=u;
elseif isa(v,'LEVELSET') && isa(u,'double')
ls=v;
v=u;
else
    error('multiplication pas possible')
end

if numel(v)>1
    error('multiplication pas possible')
end

if iseval(ls)
ls.value = v*ls.value;
else
ls.sign = v*ls.sign;
end