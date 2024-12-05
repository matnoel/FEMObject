function [B,result] = multisvd_inverse_operator(A,varargin)
% function [B,result] = multisvd_inverse_operator(A,varargin)

dim=A.dim;

solver = SEPSOLVER(dim,varargin{:});
param = getparam(solver);

if ischarin('C',varargin)
    C=getcharin('C',varargin);
else
    C=A*A';
end

if ischarin('symmetry',varargin)
    symdim = getcharin('symmetry',varargin);
    sym = zeros(1,A.dim);
    sym(symdim) = 1;
else
    sym = zeros(1,A.dim);
end
sym = logical(sym);

if ischarin('diagonal',varargin)
    diadim = getcharin('diagonal',varargin);
    dia = zeros(1,A.dim);
    dia(diadim) = 1;
else
    dia = zeros(1,A.dim);
end
dia = logical(dia);

if ischarin('initialguess',varargin)
    B=getcharin('initialguess',varargin);
else
    B=SEPMATRIX(A.dim);
end

if ischarin('W0',varargin)
    W0=getcharin('W0',varargin);
end

I=sepopereye(size(A));
switch param.errorindicator
    case {'residual','reference'}
        normref = norm(I);
end

erriter=zeros(1,param.maxiter);
errorder=zeros(1,param.maxorder);

for zz=1:param.maxorder
    if param.display
        fprintf('order #%d \n',zz)
    end
    if ischarin('W0',varargin)
        W=W0;
    else
        % W=sepoperrand(size(A));
        W=I;
    end
    alpha=W.alpha;
    for yy = 1:param.maxiter
        alpha0=alpha;
        for k=1:dim
            if B.m>0
                AWk=partialprodscal(A'-B*C,W,k);
            else
                AWk=partialprodscal(A',W,k);
            end
            AWk=reduce(AWk,k);
            
            Pk=partialprodscal(W*C,W,k);
            Pk.F(:,k)=C.F(:,k);
            Pk=reduce(Pk,k);
            
            if sym(k)
                W.F{1,k}=mylyap(Pk',Pk,-(AWk+AWk'));
            elseif dia(k)
                szf = size(AWk);
                lambda = zeros(szf(2),1);
                for ii=1:szf(2)
                    AWkii = AWk(ii,:);
                    lambda(ii) = (AWkii*Pk(:,ii)) / (AWkii*AWkii');
                end
                W.F{1,k} = spdiags(lambda,0,szf(2),szf(2));
            else
                W.F{1,k} = AWk/Pk;
            end
            
            W=normalizefuns(W);
        end
        alpha=W.alpha;
        erriter(yy)=abs(alpha-alpha0)/(alpha+alpha0);
        if param.display
            fprintf('  iteration #%d - stagnation = %d\n',yy,erriter(yy))
        end
        if (yy>1) && (erriter(yy)<param.itercrit)
            break
        end
    end
    B=B+W;
    
    switch param.errorindicator
        case {'residual','reference'}
            errorder(zz)=norm(I-B*A)/normref;
        otherwise
            errorder(zz)= sqrt(min(B.alpha)^2/sum(B.alpha.^2));
    end
    if param.display
        fprintf('  order #%d - error = %d \n',zz,errorder(zz))
    end
    
    if errorder(zz)<param.tol
        break
    end
end

result.error = errorder;

end
