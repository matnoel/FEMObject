function w = times(u,v)


if isa(u,'TIMEMATRIX') && ~isa(v,'TIMEMATRIX')
if ~isa(u.value,'cell')
   if all(size(v)==1)
    w=u;
    w.value = u.value*v;
   elseif all(u.s==1) && all(size(v)==1)
    w=u;
    w.value = v(:)*u.value;
    w.s = size(v);   
   elseif all(u.s==1)
    w = PCRADIALMATRIX(v,size(v),u); 
    return
   else
    w=u;
    w.value = times(u.value,repmat(v(:),1,length(u)));
  
   end
   if israndom(w.value)
     w = PCTIMEMATRIX(w.value,u.TIMEMODEL,w.s);  
   end
else
    error('pas programme')
end

 
elseif isa(v,'TIMEMATRIX') && ~isa(u,'TIMEMATRIX')
if ~isa(v.value,'cell')
  if all(size(u)==1)
    w=v;
    w.value = u*v.value;
  elseif all(v.s==1)  &&  all(size(u)==1)
    w=v;
    w.value = u(:)*v.value;
    w.s = size(u);  
  elseif all(v.s==1) 
    w = TIMERADIALMATRIX(u,size(u),v);
    return
  else
    w=v;
    w.value = times(repmat(u(:),1,length(v)),v.value); 
 
  end
   if israndom(w.value)
     w = PCTIMEMATRIX(w.value,v.TIMEMODEL,w.s);  
   end
else
    error('pas programme')
end    

elseif isa(v,'TIMEMATRIX') && isa(u,'TIMEMATRIX')
  n = length(u.TIMEMODEL);
      if isa(u.value,'MULTIMATRIX') && isa(v.value,'MULTIMATRIX') & ~all(sizem(u.value)==sizem(v.value))
          error('pas defini')
      elseif isa(u.value,'MULTIMATRIX')
          sm = sizem(u.value);
      else
          sm = sizem(v.value);
      end
      U = MULTIMATRIX(double(u.value),u.s,[n,prod(sizem(u.value))]);    
      V = MULTIMATRIX(double(v.value),v.s,[n,prod(sizem(v.value))]);    
      
      w = u;
      w.value = times(U,V);
      w.s = size(w.value);
      w.value = double(w.value);
      if isa(u.value,'MULTIMATRIX') || isa(v.value,'MULTIMATRIX')
      w.value = MULTIMATRIX(w.value,[prod(w.s),n],sm);
      end

end

