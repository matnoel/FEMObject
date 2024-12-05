function u = subsref(u,s)
% function u = subsref(u,s)

v = subsref(u.MULTIMATRIX,s);
if isa(v,'MULTIMATRIX') && strcmp(s.type,'()')
    u.MULTIMATRIX = v;
    
    if ~isempty(u.ximasse) && ~all(prod(sizem(u.ximasse))==prod(size(u)));
        u.ximasse = {};
    end
    
    if ~isempty(u.ximasse)
        if length(s)==1
            u.ximasse = getmultimatrix(u.ximasse,s.subs{1});
        else
            error('pas programme : ximasse pour deux subsref')
        end
    end
    
else
    u = v;
end


