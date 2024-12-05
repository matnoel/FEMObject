function ok = israndom(ls)
% function ok = israndom(ls)
% ok = 1 si la LEVELSET est aleatoire
% ok = 0 sinon
%
% See also LEVELSETS/israndom

if iseval(ls)
    ok = israndom(ls.value);
else
    if isa(ls.value{1,1},'function_handle')
        levels = ls.value{1,2};
        ok = 0 ;
        for k=1:length(levels)
            ok = ok | israndom(levels{k});
        end
    else
        ok = israndomparam(ls.value);
    end
end


