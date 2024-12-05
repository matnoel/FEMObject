function u=enrichvector(BC,u,k,choix)
rep = cell(ndims(u),1);
rep(:)={':'};

if nargin<3 | isempty(k)
    k=1;
end

if nargin<4
rep(k)={BC.ddlfree};
elseif strcmp(choix,'free')
rep(k)={getddlenrich(BC,'free')};
elseif strcmp(choix,'bloque')
rep(k)={getddlenrich(BC,'bloque')};
end

u = u(rep{:});
