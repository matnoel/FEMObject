function v = subsref(rad,s)
% function v = subsref(rad,s)

if length(s)==1 && (strcmp(s.type,'()'))
    v = rad;
    v.V = subsref(v.V,s);
elseif length(s)==1 && (strcmp(s.type,'{}'))
    
    if length(s.subs)==1 && length(s.subs{1})==1
        v = getmatrixatstep(rad,s.subs{1});
    else
        error('pas defini')
    end
end