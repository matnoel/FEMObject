function w = horzcat(u,v,varargin)
% function w = horzcat(u,v,varargin)

if nargin==1
    w = u;
else
    if isa(u,'TIMERADIALMATRIX') && isa(v,'TIMERADIALMATRIX')
        
        if u.m==v.m && norm(u.L-v.L)<eps;
            w = u;
            w.V = horzcat(u.V,v.V);
        else
            vz = sparse(size(v,1),size(v,2));
            uz = sparse(size(u,1),size(u,2));
            u.V = horzcat(u.V,vz);
            v.V = horzcat(uz,v.V);
            w = u;
            w.V = multihorzcat(u.V,v.V);
            w.L = vertcat(u.L,v.L);
            w.m = numel(w.L);
            w.D = speye(w.m);
        end
    elseif isa(u,'TIMEMATRIX')
        w = horzcat(expand(u),v);
    elseif isa(v,'TIMEMATRIX')
        w = horzcat(u,expand(v));
    elseif isa(u,'double')
        if normest(u)==0
            w = v;
            w.V = horzcat(u,w.V);
        else
            u = u * one(gettimemodel(v));
            w = horzcat(u,v);
        end
        
    elseif isa(v,'double')
        if normest(v)==0
            w = u;
            
            w.V = horzcat(w.V,v);
        else
            v = v * one(gettimemodel(u));
            w = horzcat(u,v);
        end
    else
        error('pas d�fini');
    end
    
    if nargin>2
        w = horzcat(w,varargin{:});
    end
    
    
end