function u=bloquevector(BC,u,k)
rep = cell(ndims(u),1);
rep(:)={':'};

if nargin<3
    k=1;
end

rep(k)={BC.ddlbloque};

u = u(rep{:});
