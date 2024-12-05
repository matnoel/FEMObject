function w=vertcat(u,v,varargin)


if nargin==1
    w=u;
else
if ~isa(u,'PCRADIALMATRIX') && ~isa(v,'PCRADIALMATRIX')
   w = vertcat(u,v);  
    
elseif isa(u,'PCRADIALMATRIX') && isa(v,'PCRADIALMATRIX')
  if u.m==v.m && norm(u.L-v.L)<eps;
    w=u;
    w.V = vertcat(w.V,v.V); 
 else
    vz = sparse(size(v,1),size(v,2));
    uz = sparse(size(u,1),size(u,2));
    u.V = vertcat(u.V,vz);
    v.V = vertcat(uz,v.V);
    w=u;
    
    w.V = multivertcat(u.V,v.V); %%%%%%% changement
    w.L = vertcat(u.L,v.L);
    w.m = numel(w.L);
    w.D = speye(w.m);
 end      
elseif isa(u,'PCMATRIX')
 w = vertcat(u,expand(v));
elseif isa(v,'PCMATRIX')
 w = vertcat(expand(u),v);
elseif isa(u,'double')
 if normest(u)==0
    w=v;
    w.V = vertcat(u,w.V);
 else
    u = u * one(getPC(v));  
    w = vertcat(u,v);
 end
 
elseif isa(v,'double')
 if normest(v)==0
    w=u;
    w.V = vertcat(w.V,v); 
 else
    v = v * one(getPC(u));
    w = vertcat(u,v);
 end
else
  w = vertcat(u,v); 
end

if nargin>2
  w = vertcat(w,varargin{:});
else
     if isa(w,'PCRADIALMATRIX')
   w.V = reshapem(w.V,[getm(w),1]);
   w.L = reshape(w.L,[getm(w),1]);
    end

end

end

