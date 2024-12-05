function u = subsasgn(u,s,v)
% function u = subsasgn(u,s,v)

if length(s)==1 && (strcmp(s.type,'{}'))
    
    if length(s.subs)~=1
        error('un seul indice')
    end
    
    if s.subs{1}==0
        u = setphi0(u,v);
    else
        u = setphi(u,v,s.subs{1});
    end
    
else
    error('subsref pas defini')
end
