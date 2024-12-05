function w = mtimes(u,v)

if isa(u,'PCTIMEMATRIX') & ~istime(v)
if ~isa(u.value,'cell')
   if all(size(v)==1)
    w=u;
    w.value = u.value*v;
   elseif all(u.s==1)
    w=u;
    w.value = v(:)*u.value;
    w.s = size(v);   
   elseif size(u,1)==1
    w=u;
    w.value = transpose(tranpose(u.value)*v);
    w.s = [size(u,1),size(v,2)]; 
   else
    error('pas programme')  
   end
 
else
w=u;    
for i=1:length(u.value)
    w.value{i} = u.value{i}*v;
end
w.s = size(w.value{1}); 

end

 
elseif isa(v,'PCTIMEMATRIX') & ~istime(u)
    
if ~isa(v.value,'cell')
  if all(size(u)==1)
    w=v;
    w.value = u*v.value;
  elseif all(v.s==1) 
    w=v;
    w.value = u(:)*v.value;
    w.s = size(u);   
  elseif size(v,2)==1
    w=v;
    w.value = u*v.value;
    w.s = [size(u,1),size(v,2)];     
  else
      fjdf
    error('pas programme')  
   
 
  end

else
w=v;    
for i=1:length(v.value)
    w.value{i} = u*v.value{i};
end
w.s = size(w.value{1}); 

end

elseif isa(v,'PCTIMEMATRIX') & isa(u,'PCTIMEMATRIX')
    error('pas programme')
else
    error('pas programme')
end

