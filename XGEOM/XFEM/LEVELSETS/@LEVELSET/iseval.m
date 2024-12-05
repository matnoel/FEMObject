function rep = iseval(ls)
% function rep = iseval(ls)

rep = strcmpi(class(ls),'LEVELSET') ;

if rep && isa(ls.value,'cell') && isa(ls.value{1,1},'function_handle')
    rep = 0;
end