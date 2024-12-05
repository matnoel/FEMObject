function u = subsasgn(u,s,v)
% function u = subsasgn(u,s,v)

if length(s.subs)>1
    error('only one subindex')
end
switch s.type
    case '{}'
        if isa(s.subs{1},'double') && length(s.subs{1})==1
            if isa(v,'LEVELSET') || isa(v,'LSCRACK')
                u.LS{s.subs{1}} = v;
            else
                error('right argument must be a LEVELSET')
            end
        else
            error('subindex must be of length 1')
        end
    case '()'
        if length(s.subs)==1
            if isa(s.subs{1},'double')
                if isa(v,'cell')
                    u.LS(s.subs{1}) = v;
                elseif isa(v,'LEVELSET') || isa(v,'LSCRACK')
                    u.LS(s.subs{1}) = {v};
                elseif isa(v,'LEVELSETS')
                    u.LS(s.subs{1}) = v.LS;
                else
                    error('right argument must be a LEVELSETS, a LEVELSET or a cell')
                end
            end
        end
end

u.n = length(u.LS);

u = setnumber(u);
