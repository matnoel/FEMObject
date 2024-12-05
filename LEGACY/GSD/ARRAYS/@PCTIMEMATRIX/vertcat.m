function w=vertcat(u,v,varargin)
if nargin==1
    w=u;
else
if isa(u,'PCTIMEMATRIX') & isa(v,'PCTIMEMATRIX')

 if u.s(2)>1 || v.s(2)>1
     error('pas programme');
 end
 w=u;
 w.value = vertcat(u.value,v.value); 
 w.s(1) = size(w.value,1); 

elseif ~istime(u) & isa(v,'PCTIMEMATRIX')
  u = TIMEMATRIX(u(:)*one(v.TIMEMODEL),v.TIMEMODEL,size(u));
  w = vertcat(u,v);
elseif ~istime(v) & isa(u,'PCTIMEMATRIX')  
  v = TIMEMATRIX(v(:)*one(u.TIMEMODEL),u.TIMEMODEL,size(v));
  w = vertcat(u,v); 
elseif ~israndom(u) & isa(v,'PCTIMEMATRIX')
  w = vertcat(u*one(getPC(v)),v);
elseif ~israndom(v) & isa(u,'PCTIMEMATRIX')
  w = vertcat(u,v*one(getPC(u))); 
elseif ~isa(u,'PCTIMEMATRIX') & ~isa(v,'PCTIMEMATRIX')
  w = vertcat(u,v); 
else
error('pas défini');    
end

if nargin>2
  w = vertcat(w,varargin{:});
end

end