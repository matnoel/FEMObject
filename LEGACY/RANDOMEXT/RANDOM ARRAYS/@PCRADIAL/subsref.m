function u = subsref(rad,s)
% function u = subsref(rad,s)

if length(s)==1 && (strcmp(s.type,'()'))
    if length(s.subs)==1
        u = rad.L{s.subs{1}};
    end
elseif length(s)==1 && (strcmp(s.type,'{}'))
    
    if length(s.subs)==1
        u = rad.V{s.subs{1}};
    end
    
end