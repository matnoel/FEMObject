function v=getcell(u)
if iscell(u)
if isa(u.value,'cell')
v=u.value;
elseif isa(u.value,'MULTIMATRIX') 
    v = mat2cell(u.value)
    v = getvalue(v);
else
    error('pas prevu')
end
else
    error('')
end