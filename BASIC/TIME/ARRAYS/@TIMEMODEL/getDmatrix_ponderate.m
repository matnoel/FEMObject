function Dt = getDmatrix_ponderate(L,a,varargin)
% function Dt = getDmatrix(L,a,varargin)
% Dt : matrice de deriviation numerique avec ponderation
% int(v(t)a(t)d(u(t))/dt) = int(v(t)a(t)w(t)) forall v
% -> en algebrique : w = Dt*u
%
% varargin : options
% 'ini' : prise en compte de la condition initiale, matrice augmentee
% 'basic' : int(v(t)a(t)d(u(t))/dt) = v^T Dt u en algebrique
%

dt = getdt(L);
nt = getnt(L);

switch getapproxparam(L,'type')
    case {'default','CG'}
        error('pas defini')
        switch  getapproxparam(L,'p')
            
            otherwise
                
                
        end
        
    case 'DG'
        
        switch  getapproxparam(L,'p')
            
            case 0
                a = getvalue(a);
                a = a(:);
                Dt = spdiags(a(:),0,nt,nt) + spdiags(-a(2:end),-1,nt,nt);
                
            case {1,2}
                error('pas defini')
                
                
            otherwise
                error('pas defini')
                
        end
        
        if ischarin('ini',varargin)
            u0 = sparse(size(Dt,1),1);
            u0(1)=-a(1);
            Dt = [u0,Dt];
        end
        
        
        
    otherwise
        error('pas defini')
        
        
end


if ~ischarin('basic',varargin)
    Mt = calc_masse(a);
    Dt = Mt\Dt;
end

