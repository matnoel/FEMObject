function u = subsasgn(u,s,v)
% function u = subsasgn(u,s,v)

if length(s.subs)>1
    error('only one subindex')
end
switch s.type
    case '{}'
        if isa(s.subs{1},'double') && length(s.subs{1})==1
            if isa(v,'RANDVAR')
                u.RV{s.subs{1}} = v;
            else
                error('right argument must be a RANDVAR')
            end
        else
            error('subindex must be of length 1')
        end
    case '()'
        if length(s.subs)==1
            if isa(s.subs{1},'double')
                if isa(v,'cell')
                    u.RV(s.subs{1}) = v;
                elseif isa(v,'RANDVAR')
                    u.RV(s.subs{1}) = {v};
                elseif isa(v,'RANDVARS')
                    u.RV(s.subs{1}) = v.RV;
                else
                    error('right argument must be a RANDVARS, a RANDVAR or a cell')
                end
            end
        end
end

u.M = length(u.RV);

u = setnumber(u);
u = unique(u);