function w = mrdivide(u,v)

if isa(u,'TIMEMATRIX') & ~isa(v,'TIMEMATRIX')
if ~isa(u.value,'cell')
   if all(size(v)==1)
    w=u;
    w.value = mrdivide(u.value,v);
   else
    error('pas programme')  
   end
 
else
    error('pas programme')
end

 
elseif isa(v,'TIMEMATRIX') & ~isa(u,'TIMEMATRIX')
    error('pas programme')  
    
elseif isa(v,'TIMEMATRIX') & isa(u,'TIMEMATRIX')
    error('pas programme')
    
end

