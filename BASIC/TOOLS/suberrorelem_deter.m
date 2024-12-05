function [errorelem, refelem] = suberrorelem_deter(elem,xnode,subgauss,ue,ueref,detJ,varargin)

errorelem = 0;
refelem = 0;

if ischarin('cont',varargin)
    for i=1:length(subgauss.w)
        errorelem = errorelem + error_cont(subgauss.coord(i,:),elem,xnode,ue,ueref)*abs(detJ)*subgauss.w(i);
        refelem = refelem + error_cont_ref(subgauss.coord(i,:),elem,xnode,ueref)*abs(detJ)*subgauss.w(i);
    end
elseif ischarin('depla',varargin) || isempty(varargin)
    for i=1:length(subgauss.w) 
        errorelem = errorelem + error_depla(subgauss.coord(i,:),elem,xnode,ue,ueref)*abs(detJ)*subgauss.w(i);
        refelem = refelem + error_depla_ref(subgauss.coord(i,:),elem,xnode,ueref)*abs(detJ)*subgauss.w(i);
    end
end

function rep = error_depla(gauss,elem,xnode,ue,uref)

Ne = calc_N(elem,xnode,gauss);
uapprox = Ne*ue;
u_ref = Ne*uref;
rep = double((uapprox-u_ref)'*(uapprox-u_ref));

return

function rep = error_depla_ref(gauss,elem,xnode,uref)

Ne = calc_N(elem,xnode,gauss);
u_ref = Ne*uref;
rep = double(u_ref'*u_ref);

return

function rep = error_cont(gauss,elem,xnode,ue,uref)

mat = getmaterial(elem);
se = sigma(mat,elem,xnode,gauss,ue);
se_ref = sigma(mat,elem,xnode,gauss,uref);
rep = double((se-se_ref)'*(se-se_ref));

return

function rep = error_cont_ref(gauss,elem,xnode,uref)

mat = getmaterial(elem);
se_ref = sigma(mat,elem,xnode,gauss,uref);
rep = double(se_ref'*se_ref);

return
