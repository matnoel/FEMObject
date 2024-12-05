function u = subsref(u,s)
% function u = subsref(u,s)

if length(s)==1 && (strcmp(s.type,'{}'))
    
    if length(s.subs)~=1
        error('un seul indice')
    end
    
    if s.subs{1}==0
        u = u.phi0;
    else
        u = u.phi{s.subs{1}};
    end
    
elseif length(s)==1 && (strcmp(s.type,'()'))
    
    u.phi0 = subsref(u.phi0,s);
else
    error('subsref pas defini')
end
