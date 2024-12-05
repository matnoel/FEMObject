function w=vertcat(u,v,varargin)
if nargin==1
    w=u;
else
if isa(u,'TIMERADIALMATRIX') && isa(v,'TIMERADIALMATRIX')
  if u.m==v.m && norm(u.L-v.L)<eps;
    w=u;
    w.V = vertcat(w.V,v.V); 
 else
    vz = sparse(size(v,1),size(v,2));
    uz = sparse(size(u,1),size(u,2));
    u.V = vertcat(u.V,vz);
    v.V = vertcat(uz,v.V);
    w=u;
    w.V = multihorzcat(u.V,v.V);
    w.L = vertcat(u.L,v.L);
    w.m = numel(w.L);
    w.D = speye(w.m);
 end      
elseif isa(u,'TIMEMATRIX')
 w = vertcat(u,expand(v));
elseif isa(v,'TIMEMATRIX')
 w = vertcat(expand(u),v);
elseif isa(u,'double')
 if normest(u)==0
    w=v;
    w.V = vertcat(u,w.V);
 else
    u = u * one(gettimemodel(v));  
    w = vertcat(u,v);
 end
 
elseif isa(v,'double')
 if normest(v)==0
    w=u;
    w.V = vertcat(w.V,v); 
 else
    v = v * one(gettimemodel(u));
    w = vertcat(u,v);
 end
else
error('pas défini');    
end

if nargin>2
  w = vertcat(w,varargin{:});
end

end


