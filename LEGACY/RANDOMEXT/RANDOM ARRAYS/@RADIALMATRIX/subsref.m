function [u,v] = subsref(rad,s)
% function [u,v] = subsref(rad,s)

if length(s)==1 && (strcmp(s.type,'()'))
    
    u = rad ;
    u.V = subsref(u.V,s);
    
elseif length(s)==1 && (strcmp(s.type,'{}'))
    
    if length(s.subs)==1
        u = rad.V{s.subs{1}};
        if nargout>1
            v = rad.L(s.subs{1});
        end
    end
    
    
end