function var = getvar(T,i)
if nargin ==1
    var = T.var;
elseif nargin==2 && isa(i,'double')
    nbcell = length(i);
    if nbcell==1
        var=T.Cvar{i};
    else
        var = cell(nbcell,1);
        for ii=1:nbcell
            var{ii}=T.Cvar{i(ii)};
        end
    end
end