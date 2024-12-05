function param = getparam(rv,paramname)
% function param = getparam(rv,paramname)

if nargin==1
    param = rv.param';
    param = struct(param{:});
else
    param = rv.param;
    if ~isa(paramname,'char') && ~(isa(paramname,'cell') && length(paramname)==1)
        error('indiquer le nom d''un seul parametre' )
    end
    
    [rep,pos] = ischarin(paramname,param(:,1));
    param = param{pos(find(rep)),2};
end

