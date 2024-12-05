function ke = masslscarac_elem(elem,xnodee,lse,varargin)

method = getcharin('method',varargin);
switch method
    case 'carac'
        
        ke=zerosND(getnbddl(elem),getnbddl(elem),getnbelem(elem));
        
        %gauss = calc_gauss(elem,10);
        
        if all(lse<=0)
            gauss = calc_gauss(elem,6);
        elseif ~all(lse>=0)
            [gauss,gaussout] = calc_lssubgauss(elem,lse,6);
            
        end
        
        x = calc_x(elem,xnodee,gauss.coord);
        
        N = calc_N(elem,xnodee,gauss.coord);
        DN = calc_DN(elem,xnodee,gauss.coord);
        
        detJ=calc_detJ(elem,xnodee,gauss.coord);
        
        psi = -N*lse;
        Dpsi = -DN*lse;
        
        %     psi = (N*lse).^2;
        %     Dpsi = 2*(N*lse).*(DN*lse);
        %
        %   psi =  exp(-N*lse)-1;
        %   Dpsi = -exp(-N*lse).*(DN*lse);
        %
        ke(:,:)=sum(gauss.w*(psi*DN')*(psi*DN)*abs(detJ),4)+ ...
            sum(gauss.w*(N'*(Dpsi)')*(psi*DN)*abs(detJ),4)+...
            sum(gauss.w*(psi*DN')*((Dpsi)*N)*abs(detJ),4)...
            +sum(gauss.w*(N'*(Dpsi)')*((Dpsi)*N)*abs(detJ),4);
        
        
        
end


