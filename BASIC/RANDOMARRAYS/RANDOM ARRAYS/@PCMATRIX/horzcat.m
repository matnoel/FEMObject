function w = horzcat(u,v,varargin)
% function w = horzcat(u,v,varargin)

if nargin==1
    w = u;
else
    if isa(u,'PCMATRIX') && isa(v,'PCMATRIX')
        w = u;
        w.MULTIMATRIX = horzcat(u.MULTIMATRIX,v.MULTIMATRIX);
    elseif isa(u,'double')
        if normest(u)==0
            w = v;
            w.MULTIMATRIX = horzcat(u,w.MULTIMATRIX);
        else
            u = expand(u * one(getPC(v)));
            w = horzcat(u,v);
        end
        
    elseif isa(v,'double')
        if normest(v)==0
            w = u;
            w.MULTIMATRIX = horzcat(w.MULTIMATRIX,v);
        else
            v = expand(v * one(getPC(u)));
            w = horzcat(u,v);
        end
    else
        error('pas d�fini');
    end
    
    if nargin>2
        w = horzcat(w,varargin{:});
    end
    
end