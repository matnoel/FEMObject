function u = subsref(u,s)
% function u = subsref(u,s)

switch s.type
    case '{}'
        if isa(u.value,'cell')
            u = subsref(u.value,s);
        elseif isa(u.value,'double') || israndom(u.value) || isa(u.value,'MULTIMATRIX')
            u = reshape(u.value(:,s.subs{1}),u.s);
        elseif isa(u.value,'FEELEMFIELD')
            u = u.value(:,s.subs{1});
        else
            error('pas programme')
        end
    case '()'
        if isa(u.value,'double')
            v = MULTIMATRIX(u.value,u.s,[length(gettapprox(u)),1]);
            v = subsref(v,s);
            u.s = size(v);
            u.value = double(v);
        elseif strcmp(class(u.value),'MULTIMATRIX')
            u.value = subsref(u.value,s);
            u.s = size(u.value);
        elseif isa(u.value,'FEELEMFIELD')
            
            if length(s.subs)==1
                s.subs=[s.subs {':'}];
                u.value = subsref(u.value,s);
                u.s(1) = size(u.value,1);
            else
                warning('ambiguite sur le subsref')
                u = subsref(u.value,s);
            end
            
        elseif israndom(u.value)
            if length(s.subs)==1 && u.s(2)==1
                s.subs{2} = ':';
                u.value = subsref(u.value,s);
                u.s(1) = size(u.value,1);
            else
                error('pas programme')
            end
        end
end
