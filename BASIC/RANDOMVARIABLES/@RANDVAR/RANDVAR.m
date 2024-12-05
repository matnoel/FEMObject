function rv = RANDVAR(type,param,domain)
% function rv = RANDVAR(type,param,domain)

if nargin==0
    rv.type = '';
    rv.param = cell(0,1);
    rv.number = [];
    rv = class(rv,'RANDVAR');        
    inferiorto('CONDRANDVAR')

elseif isa(type,'RANDVAR')
    rv = type;
else
    rv.type = type;
    paramnames = fieldnames(param);
    rv.param = struct2cell(param);

    rv.param = [paramnames(:),rv.param(:)];

    rv.number = [];

    rv = class(rv,'RANDVAR');        
    inferiorto('CONDRANDVAR')
end
