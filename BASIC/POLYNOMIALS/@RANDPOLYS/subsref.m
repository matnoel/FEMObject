function varargout = subsref(u,s)
% function varargout = subsref(u,s)

if length(s.subs)>1
    error('only one subindex')
end

switch s.type
    case '{}'
        
        if isa(s.subs{1},'double')
            varargout = cell(1,length(s.subs{1}));
            varargout(:) = u.h(s.subs{1});
        elseif strcmp(s.subs{1},':')
            varargout = u.h ;
        else
            error(' ')
        end
        
    case '()'
        
        u.h = u.h(s.subs{1});
        u.M = length(u.h);
        varargout{1} = u;
        
end
