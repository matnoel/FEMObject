function varargout = subsref(u,s)
% function varargout = subsref(u,s)

if length(s.subs)>1
    error('only one subindex')
end

switch s.type
    case '{}'
        
        if isa(s.subs{1},'double')
            varargout = cell(1,length(s.subs{1}));
            varargout(:) = u.RV(s.subs{1});
        elseif strcmp(s.subs{1},':')
            varargout = u.RV ;
        else
            error(' ')
        end
        
    case '()'
        
        u.RV = u.RV(s.subs{1});
        u.M = length(u.RV);
        varargout{1} = u;
        
end

