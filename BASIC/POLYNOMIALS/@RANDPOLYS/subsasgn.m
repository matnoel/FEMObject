function u = subsasgn(u,s,v)
% function u = subsasgn(u,s,v)

if length(s.subs)>1
    error('only one subindex')
end
switch s.type
    case '{}'
        if isa(s.subs{1},'double') && length(s.subs{1})==1
            if isa(v,'RANDPOLY')
                u.h{s.subs{1}} = v;
                u.M = length(u.h);
            else
                error('right argument must be a RANDPOLY')
            end
        else
            error('subindex must be of length 1')
        end
    case '()'
        if length(s.subs)==1
            if isa(s.subs{1},'double')
                if isa(v,'cell')
                    u.h(s.subs{1}) = v;
                elseif  isa(v,'RANDPOLY')
                    u.h(s.subs{1}) = {v};
                elseif isa(v,'RANDPOLYS')
                    u.h(s.subs{1}) = v.h;
                else
                    error('right argument must be a RANDPOLYS, a RANDPOLY or a cell')
                end
            end
        end
end

u.M = length(u.h);

u = setnumber(u);