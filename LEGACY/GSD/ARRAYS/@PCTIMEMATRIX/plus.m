function w = plus(u,v)

if isa(u,'PCTIMEMATRIX') & ~isa(v,'PCTIMEMATRIX')
w=u;
if isa(u.value,'cell') 
if  isa(v,'TIMEMATRIX')
for i=1:length(u.value)
u.value{i} = u.value{i} + v{i} ;
end  
else
for i=1:length(u.value)
u.value{i} = u.value{i} + v ;
end
u.s = size(u.value{1});
end
    
else
if isa(v,'TIMEMATRIX')
  w = plus(u,v*one(getPC(u.value)));
elseif isa(v,'double')   
if all(size(v)==1)
  w=u;
  w.value = w.value + one(getPC(u.value));
else 
  v = TIMEMATRIX(repmat(v(:),1,length(u.TIMEMODEL)),u.TIMEMODEL,size(v));
  w = plus(u,v);
end
else    
    error('pas programme')
end
    
end

elseif isa(v,'PCTIMEMATRIX') & ~isa(u,'PCTIMEMATRIX')

 w = plus(v,u);   
    
elseif isa(v,'PCTIMEMATRIX') & isa(u,'PCTIMEMATRIX')

    
if isa(u.value,'cell') && isa(v.value,'cell')
w=u;
for i=1:length(u.value)
w.value{i} = u.value{i} + v.value{i} ;
end  
w.s = size(w.value{1});
elseif ~isa(u.value,'cell') && ~isa(v.value,'cell')
w=u;
if all(size(u)==size(v))
    w.value = u.value+v.value;
else
error('pas programme');
end

elseif isa(u.value,'cell')
w = plus(u,mat2cell(v));
elseif isa(v.value,'cell')
w = plus(mat2cell(u),v);
end


end

