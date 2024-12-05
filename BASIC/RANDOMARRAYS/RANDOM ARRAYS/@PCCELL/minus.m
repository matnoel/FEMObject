function w=minus(u,v)

if isa(v,'PCCELL') & isa(u,'PCCELL')
    w=u;
    s=size(w.value{1});
    u=getmatrix(u.value);
    v=getmatrix(v.value);
    w.value = setmatrix(u-v,s(1),s(2));

else
       error('minus not defined')     
end
    
