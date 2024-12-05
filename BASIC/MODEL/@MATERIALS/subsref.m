function varargout = subsref(u,s)
% function varargout = subsref(u,s)

if length(s.subs)>1
    error('only one subindex')
end

switch s.type
    case '{}'
        
        if isa(s.subs{1},'double')
            varargout = cell(1,length(s.subs{1}));
            varargout(:) = u.MAT(s.subs{1});
        elseif strcmp(s.subs{1},':')
            varargout = u.MAT;
        else
            error(' ')
        end
        
    case '()'
        
        u.MAT = u.MAT(s.subs{1});
        u.n = length(u.MAT);
        varargout{1} = u;
        
end
