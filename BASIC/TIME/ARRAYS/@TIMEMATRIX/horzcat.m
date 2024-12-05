function w = horzcat(u,v,varargin)
% function w = horzcat(u,v,varargin)

if nargin==1
    w = u;
else
    
    if isa(u,'TIMEMATRIX') && isa(v,'TIMEMATRIX')
        
        U = MULTIMATRIX(u.value,u.s,[length(u.TIMEMODEL),1]);
        V = MULTIMATRIX(v.value,v.s,[length(u.TIMEMODEL),1]);
        W = horzcat(U,V);
        w = u;
        w.s = size(W);
        w.W = double(W);
        
        
    elseif ~isa(u,'TIMEMATRIX') && isa(v,'TIMEMATRIX')
        
        u = TIMEMATRIX(u*one(v.TIMEMODEL),v.TIMEMODEL,size(u));
        w = horzcat(u,v);
        
    elseif isa(u,'TIMEMATRIX') && ~isa(v,'TIMEMATRIX')
        v = TIMEMATRIX(v*one(u.TIMEMODEL),u.TIMEMODEL,size(v));
        w = horzcat(u,v);
        
    else
        error('pas defini');
    end
    
    if nargin>2
        w = horzcat(w,varargin{:});
    end
    
end