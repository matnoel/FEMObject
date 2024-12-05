function varargout = subsref(u,s)
% function varargout = subsref(u,s)

if length(s.subs)>1
    error('only one subindex')
end

switch s.type
    case '{}'
        
        if isa(s.subs{1},'double')
            varargout = cell(1,length(s.subs{1}));
            varargout(:) = u.LS(s.subs{1});
        elseif strcmp(s.subs{1},':')
            varargout = u.LS ;
        else
            error(' ')
        end
        
    case '()'
        
        for i=1:u.n
            %if isalevelset(u.LS{i})
            if iseval(u.LS{i})
                u.LS{i} = subsref(u.LS{i},s);
            end
        end
        varargout = u.LS;
        
end
