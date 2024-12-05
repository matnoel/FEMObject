function w = times(u,v)


if isa(u,'PCTIMEMATRIX') & ~istime(u)
if ~isa(u.value,'cell')
   if all(size(v)==1)
    w=u;
    w.value = u.value*v;
   elseif all(u.s==1)
    w=u;
    w.value = v(:)*u.value;
    w.s = size(v);   
   else
    w=u;
    w.value = times(u.value,repmat(v(:),1,length(u.TIMEMODEL)));
  
   end
 
else
    error('pas programme')
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
  else
    w=v;
    w.value = times(repmat(u(:),1,length(v.TIMEMODEL)),v.value); 
    
  end

else
    error('pas programme')
end    

elseif isa(v,'PCTIMEMATRIX') & isa(u,'PCTIMEMATRIX')
  w=u;
  w.value = u.value.*v.value;
  if all(size(u)==1)
  w.s = v.s;
  elseif all(size(v)==1)
   w.s = u.s;
  end
else
error('pas programme')
end

