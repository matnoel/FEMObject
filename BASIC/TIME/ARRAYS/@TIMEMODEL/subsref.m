function a = subsref(T,s)
% function a = subsref(T,s)

if length(s)>1
    error('subsref pas defini')
end

switch s.type
    case '.'
        a = getfield(struct(T),s.subs);
    case '()'
        if isa(s.subs{1},'function_handle')
            t = gettapprox(T);
            a = TIMEMATRIX(double(s.subs{1}(t)),T);
        else
            a = T.t(s.subs{:});
        end
        
    otherwise
        error('subsref pas defini')
end


