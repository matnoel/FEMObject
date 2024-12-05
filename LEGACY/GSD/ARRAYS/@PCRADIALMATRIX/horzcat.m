function w = horzcat(u,v,varargin)
% function w = horzcat(u,v,varargin)

if nargin==1
    w = u;
else
    if ~isa(u,'PCRADIALMATRIX') && ~isa(v,'PCRADIALMATRIX')
        w = horzcat(u,v);
    elseif   isa(u,'PCRADIALMATRIX') && isa(v,'PCRADIALMATRIX')
        
        if u.m==v.m && norm(u.L-v.L)<eps;
            w = u;
            w.V = horzcat(u.V,v.V);
        else
            vz = sparse(size(v,1),size(v,2));
            uz = sparse(size(u,1),size(u,2));
            u.V = horzcat(u.V,vz);
            v.V = horzcat(uz,v.V);
            w = u;
            try
                w.V = multihorzcat(u.V,v.V);
            catch
                w.V = multivertcat(u.V,v.V);
            end
            w.L = vertcat(u.L,v.L);
            w.m = numel(w.L);
            w.D = speye(w.m);
        end
    elseif isa(u,'PCMATRIX')
        w = horzcat(u,expand(v));
    elseif isa(v,'PCMATRIX')
        w = horzcat(expand(u),v);
    elseif isa(u,'double')
        if normest(u)==0
            w = v;
            w.V = horzcat(u,w.V);
        else
            u = u * one(getPC(v));
            w = horzcat(u,v);
        end
        
    elseif isa(v,'double')
        if normest(v)==0
            w = u;
            
            w.V = horzcat(w.V,v);
        else
            v = v * one(getPC(u));
            w = horzcat(u,v);
        end
    else
        w = horzcat(u,v);
    end
    
    if nargin>2
        w = horzcat(w,varargin{:});
    else
        if isa(w,'PCRADIALMATRIX')
            w.V = reshapem(w.V,[getm(w),1]);
            w.L = reshape(w.L,[getm(w),1]);
        end
        
        
    end
    
    
end
