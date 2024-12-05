function u = subsasgn(u,s,b)
% function u = subsasgn(u,s,b)

switch class(b)
    case 'POINT'
        switch s.type
            case '()'
                if length(s.subs)==1
                    u.MYDOUBLEND(:,:,s.subs{1})=b.MYDOUBLEND(:,:,:);
                end
        end
end
