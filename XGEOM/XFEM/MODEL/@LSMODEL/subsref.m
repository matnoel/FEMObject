function a = subsref(S,s)
% function a = subsref(S,s)

if strcmp(s(1).type,'.') && strcmp(s(1).subs,'ls')
    a = S.ls;
    if length(s)>1
        a = subsref(a,s(2:end));
    end
else
    a = subsref(S.MODEL,s);
end