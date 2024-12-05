function u = subsasgn(u,s,v)
% function u = subsasgn(u,s,v)

if length(s.subs)>1
    error('only one subindex')
end
switch s.type
    case '{}'
        if isa(s.subs{1},'double') && length(s.subs{1})==1
            if isa(v,'MATERIAL')
                u.MAT{s.subs{1}} = v;
            else
                error('right argument must be a MATERIAL')
            end
        else
            error('subindex must be of length 1')
        end
    case '()'
        if length(s.subs)==1
            if isa(s.subs{1},'double')
                if isa(v,'cell')
                    u.MAT(s.subs{1}) = v;
                elseif isa(v,'MATERIAL')
                    u.MAT(s.subs{1}) = {v};
                elseif isa(v,'MATERIALS')
                    u.MAT(s.subs{1}) = v.MAT;
                else
                    error('right argument must be a MATERIALS, a MATERIAL or a cell')
                end
            end
        end
end

u.n = length(u.MAT);

u = setnumber(u);
