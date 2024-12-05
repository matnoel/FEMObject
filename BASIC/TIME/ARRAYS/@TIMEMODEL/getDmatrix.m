function Dt = getDmatrix(L,varargin)
% function Dt = getDmatrix(L,varargin)
% Dt : matrice de deriviation numerique associée a un solveur
% int(v(t)d(u(t))/dt) = int(v(t)w(t)) forall v
% -> en algebrique : w = Dt*u
%
% varargin : options
% 'ini' : prise en compte de la condition initiale, matrice augmentee
% 'basic' : int(v(t)d(u(t))/dt) = v^T Dt u en algebrique
%

dt = getdt(L);
nt = getnt(L);

switch getapproxparam(L,'type')
    case {'default','CG'}
        switch  getapproxparam(L,'p')
            case 0
                Dt = spdiags(ones(nt,1),0,nt,nt) + spdiags(-ones(nt-1,1),-1,nt,nt);
            case 1
                dtd=-1/2*[ones(nt,1);0];
                dtu=1/2*[0;ones(nt,1)];
                dtm = dtd+dtu;
                
                Dt = spdiags([dtd,dtm,dtu],[-1,0,1],nt+1,nt+1);
            otherwise
                error('pas defini')
        end
    case 'DG'
        switch  getapproxparam(L,'p')
            case 0
                Dt = spdiags(ones(nt,1),0,nt,nt) + spdiags(-ones(nt-1,1),-1,nt,nt);
            case {1,2}
                p=getapproxparam(L,'p') ;
                if p==1
                    matDelem = [-1/2,1/2;-1/2,1/2];
                elseif p==2
                    matDelem = [-1/2,2/3,-1/6;-2/3,0,2/3;1/6,-2/3,1/2];
                end
                
                I=[];
                J=[];
                V=[];
                for i=1:p+1
                    for j=1:p+1
                        Iadd = [i:p+1:nt*(p+1)]';
                        Jadd = [j:p+1:nt*(p+1)]';
                        I = [I;Iadd];
                        J = [J;Jadd];
                        V = [V;ones(length(Iadd),1)*matDelem(i,j)];
                    end
                end
                Ipadd = [1:p+1:nt*(p+1)]';
                Imadd = [p+2:p+1:nt*(p+1)]';
                Jmadd = [p+1:p+1:(nt-1)*(p+1)]';
                
                I = [I;Ipadd];
                J = [J;Ipadd];
                V = [V;ones(length(Ipadd),1)];
                
                I = [I;Imadd];
                J = [J;Jmadd];
                V = [V;-ones(length(Imadd),1)];
                
                Dt = sparse(I,J,V);
            otherwise
                error('pas defini')
        end
        
        if ischarin('ini',varargin)
            u0 = sparse(size(Dt,1),1);
            u0(1)=-1;
            Dt = [u0,Dt];
        end
        
    otherwise
        error('pas defini')
end


if ~ischarin('basic',varargin)
    Mt = getMmatrix(L);
    Dt = Mt\Dt;
end

