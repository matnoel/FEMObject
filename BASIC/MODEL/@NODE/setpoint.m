function u = setpoint(u,P)

if ~strcmp(class(P),'POINT')
    error('rentrer un POINT')
end

if getnbnode(u)~=numel(P)
    error('rentrer le bon nombre de point')
end

u.POINT = P;
