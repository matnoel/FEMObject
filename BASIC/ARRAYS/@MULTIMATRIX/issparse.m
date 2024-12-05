function a=issparse(u)
if isa(u.value,'cell')
a = issparse(u.value{1});
else
a=issparse(u.value);
end