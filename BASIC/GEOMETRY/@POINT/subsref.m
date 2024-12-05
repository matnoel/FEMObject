function u = subsref(u,s)
% function u = subsref(u,s)

switch s.type
    case '.'
        u = get(u,s.subs);
    case {'()'}
        switch length(s.subs)
            case 1
                u.MYDOUBLEND = u.MYDOUBLEND(:,:,s.subs{1});
                ss = squeeze(s.subs{1});
                
                if ~isa(s.subs{1},'char') && (size(ss,2)>1 || size(ss,2)==0 || ndims(ss)>2 || ndims(ss)==0)
                    u.MYDOUBLEND = reshape(u.MYDOUBLEND,[size2D(u.MYDOUBLEND),size(s.subs{1})]);
                end
            otherwise
                u.MYDOUBLEND = u.MYDOUBLEND(:,:,s.subs{:});
        end
    otherwise
        error('non defini')
        
end

