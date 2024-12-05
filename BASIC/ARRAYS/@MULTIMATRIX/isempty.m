function x=isempty(a)
if isa(a.value,'cell')
 x= any(size(a.value)==0)| any(a.s==0);
else
    x=any(size(a.value)==0);
end