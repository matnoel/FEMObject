function [norm_error, norm_ref] = calc_errorfield(M,u,fun,varargin)
% [norm_error, norm_ref] = calc_errorfield(M,u,fun,varargin)
% varargin : 'type'
% 'type' peut �tre 'depla' ou 'cont'
% 'Mref' MODEL finalis� : utilis� dans le cas d'une solution de r�f�rence
% qui n'est pas analytique
% 'Mrefcalcul' autre MODEL finalis� pour un calcul de l'erreur sur un autre maillage (plus fin ou plus grossier) 
if isa(fun,'function_handle')
    if israndom(u)
        u = unfreevector(M,u);
        Mref = getcharin('Mref',varargin,M);
        [norms] = calc_elemfields(Mref,@calc_error,M,u,fun,varargin{:});
        norm_error=norms{1};
        norm_ref = norms{2};
    else
        u = unfreevector(M,u);
        [norms] = calc_elemfields(M,@calc_error_deter,M,u,fun,varargin{:});
        norm_error=norms{1};
        norm_ref = norms{2};    
    end
else
    Mref = getcharin('Mref',varargin);
    Mrefcalcul = getcharin('Mrefcalcul',varargin,Mref);
    if israndom(u)
        u1 = unfreevector(M,u);
        u1 = transfer(M,Mref,u1);
        u1 = project(u1,getPC(fun));
        uref = unfreevector(Mref,fun);
        [norms] = calc_elemfields(Mrefcalcul,@calc_error_Mref,Mref,uref,u1,varargin{:});
        norm_error=norms{1};
        norm_ref = norms{2};
    else
        u1 = unfreevector(M,u);
        u1 = transfer(M,Mref,u1);
        uref = unfreevector(Mref,fun);
        [norms] = calc_elemfields(Mref,@calc_error_deter,Mref,uref,u1,varargin{:});
        norm_error=norms{1};
        norm_ref = norms{2};    
    end
end    
