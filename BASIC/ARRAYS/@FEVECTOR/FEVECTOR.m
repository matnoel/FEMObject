function v = FEVECTOR(a,M,fedim)
% function v = FEVECTOR(a,M,fedim)
% a : array
% M : MODEL ou FEARRAY ou double pour indiquer les ddlbloque
% fedim : donne la dimension element fini (1 par defaut)

if nargin==1 & isa(varargin{1},'FEVECTOR')
    v=a;

else
    if nargin==1
        ddlbloque = [];
    elseif isa(M,'MODEL') | isa(M,'FEARRAY')
        ddlbloque = getddlbloque(M);
    else
        ddlbloque = M ;   
    end

    if nargin==2
        fedim=1;
    end

    a = FEARRAY(a,fedim,ddlbloque);
    v=struct();
    v=class(v,'FEVECTOR',a);
    superiorto('MYDOUBLE')
    inferiorto('PCARRAY')

end
