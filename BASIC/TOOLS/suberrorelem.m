function [errorelem refelem] = suberrorelem(elem,xnode,subgauss,ue,fun,detJ,varargin)

errorelem = 0;
refelem = 0;

if ischarin('depla',varargin) || isempty(varargin)
    if isa(fun,'MYDOUBLEND')
        for i=1:length(subgauss.w)
            errorelem = errorelem + error_depla(subgauss.coord(i,:),elem,xnode,ue,fun)*detJ*subgauss.w(i);
            refelem = refelem + error_depla_ref(subgauss.coord(i,:),elem,xnode,fun)*detJ*subgauss.w(i);   
        end
    else
        for i=1:length(subgauss.w)
            errorelem = errorelem + error_deplaex(subgauss.coord(i,:),elem,xnode,ue,fun)*detJ*subgauss.w(i);
            refelem = refelem + error_depla_refex(subgauss.coord(i,:),elem,xnode,fun)*detJ*subgauss.w(i);      
        end
    end
elseif ischarin('cont',varargin)
    if isa(fun,'MYDOUBLEND')  
        for i=1:length(subgauss.w)  
            errorelem = errorelem + error_cont(subgauss.coord(i,:),elem,xnode,ue,fun)*detJ*subgauss.w(i);
            refelem = refelem + error_cont_ref(subgauss.coord(i,:),elem,xnode,fun)*detJ*subgauss.w(i);
        end
    else
        for i=1:length(subgauss.w)
            errorelem = error_contex(subgauss.coord,elem,xnode,ue,fun)*detJ*sum(subgauss.w);
            refelem = error_cont_refex(subgauss.coord,elem,xnode,fun)*detJ*sum(subgauss.w);   
        end
    end
end

function f1 = error_depla(gauss,elem,xnode,ue,uref)

Ne = calc_N(elem,xnode,gauss);
uapprox = Ne*ue;
u_ref = Ne*uref;
f1 = double((uapprox-u_ref)'*(uapprox-u_ref));

return

function f2 = error_depla_ref(gauss,elem,xnode,uref)

Ne = calc_N(elem,xnode,gauss);
u_ref = Ne*uref;
f2 = double(u_ref'*u_ref);

return

function f1 = error_deplaex(gauss,elem,xnode,ue,fun)

[Ne detJ x] = calc_N(elem,xnode,gauss);
uapprox = Ne*ue;
uex = fun(x(1,1),x(1,2));
f1 = double((uapprox-uex)'*(uapprox-uex));

return

function f2 = error_depla_refex(gauss,elem,xnode,fun)

x = calc_x(elem,xnode,gauss);
uex = fun(x(1,1),x(1,2));
f2 = double(uex'*uex);

return

function f1 = error_cont(gauss,elem,xnode,ue,uref)

mat = getmaterial(elem);
se = sigma(mat,elem,xnode,gauss,ue);
se_ref = sigma(mat,elem,xnode,gauss,uref);
f1 = double((se-se_ref)'*(se-se_ref));

return

function f2 = error_cont_ref(gauss,elem,xnode,uref)

mat = getmaterial(elem);
se_ref = sigma(mat,elem,xnode,gauss,uref);
f2 = double(se_ref'*se_ref);

return

function f1 = error_contex(gauss,elem,xnode,ue,fun)

mat = getmaterial(elem);
se = sigma(mat,elem,xnode,gauss,ue);
x = calc_x(elem,xnode,gauss);
sigex = fun(x(1,1),x(1,2));
f1 = double((se-sigex)'*(se-sigex));

return

function f2 = error_cont_refex(gauss,elem,xnode,fun)

x = calc_x(elem,xnode,gauss);
sigex = fun(x(1,1),x(1,2));
f2 = double(sigex'*sigex);

return
