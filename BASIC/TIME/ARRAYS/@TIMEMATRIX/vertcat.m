function w=vertcat(u,v,varargin)
if nargin==1
    w=u;
else

if isa(u,'TIMEMATRIX') & isa(v,'TIMEMATRIX')

 U = MULTIMATRIX(u.value,u.s,[length(u.TIMEMODEL),1]);
 V = MULTIMATRIX(v.value,v.s,[length(u.TIMEMODEL),1]);
 W = vertcat(U,V);
 w=u;
 w.s = size(W);
 w.value = double(W);
 
elseif ~isa(u,'TIMEMATRIX') & isa(v,'TIMEMATRIX')
     
     u = TIMEMATRIX(u*one(v.TIMEMODEL),v.TIMEMODEL,size(u));
     w = vertcat(u,v);
 
elseif isa(u,'TIMEMATRIX') & ~isa(v,'TIMEMATRIX')
  
  v = TIMEMATRIX(v*one(u.TIMEMODEL),u.TIMEMODEL,size(v));
  w = vertcat(u,v);
 
else
error('pas défini');    
end

if nargin>2
  w = vertcat(w,varargin{:});
end

end