function w = plus(u,v)

if isa(u,'TIMEMATRIX') & ~isa(v,'TIMEMATRIX')
w=u;
if isa(u.value,'cell')
for i=1:length(u.value)
u.value{i} = u.value{i} + v ;
end
u.s = size(u.value{1});

else

if all(size(v)==1)
w.value = u.value + v ;
else
w.value = u.value + repmat(v(:),1,length(u.TIMEMODEL)) ;
end

    
end

elseif isa(v,'TIMEMATRIX') & ~isa(u,'TIMEMATRIX')

 w = plus(v,u);   
    
elseif isa(v,'TIMEMATRIX') & isa(u,'TIMEMATRIX')
w=u;

if all(size(u)==size(v))
    w.value = u.value+v.value;
else
    n = length(u.TIMEMODEL);
    if isa(u.value,'MULTIMATRIX') & isa(v.value,'MULTIMATRIX') & ~all(sizem(u.value)==sizem(v.value))
          error('pas defini')
    elseif isa(u.value,'MULTIMATRIX')
          sm = sizem(u.value);
    else
          sm = sizem(v.value); 
    end  
    U = MULTIMATRIX(double(u.value),u.s,[n,prod(sizem(u.value))]);    
    V = MULTIMATRIX(double(v.value),v.s,[n,prod(sizem(v.value))]);    
    w.value = U+V;
    w.s = size(w.value);
    w.value = double(w.value);
    if isa(u.value,'MULTIMATRIX') | isa(v.value,'MULTIMATRIX')
    w.value = MULTIMATRIX(w.value,[prod(w.s),n],sm);
    end

end
end

