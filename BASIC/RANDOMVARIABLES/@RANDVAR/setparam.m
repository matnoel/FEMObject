function rv = setparam(rv,paramname,paramval)
% function rv = setparam(rv,paramname,paramval)

if nargin==2 && isa(paramname,'cell')
    if ~all(size(paramname)==size(rv.param))
        error('rentrer le bon nombre de parametres')
    end
    rv.param = paramname ;
elseif nargin==3 && isa(paramname,'char')
    param = rv.param';
    param = struct(param{:});
    param = setfield(param,paramname,paramval);
    paramnames = fieldnames(param);
    param = struct2cell(param);
    rv.param = [paramnames(:),param(:)];
end