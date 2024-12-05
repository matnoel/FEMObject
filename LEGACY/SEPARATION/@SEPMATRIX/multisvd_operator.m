function [u,result] = multisvd_operator(A,varargin)
% function [u,result] = multisvd_operator(A,varargin)

A=normalizefuns(A);
dim = getdim(A);

solver = SEPSOLVER(dim,varargin{:});
param = getparam(solver);

if ischarin('w0',varargin)
    w0=getcharin('w0',varargin);
end
erriter = zeros(1,param.maxiter);
errorder = zeros(1,param.maxorder);

switch param.errorindicator
    case {'residual','reference'}
        normref = norm(A);
end

if ischarin('initialguess',varargin)
    u = getcharin('initialguess',varargin);
    A=A-u;
else
    u=SEPMATRIX(dim);
end
for zz=1:param.maxorder
    if param.display
        fprintf('order #%d \n',zz)
    end
    
    if ischarin('w0',varargin)
        w=w0;
    else
        w=sepoperrand(size(A));
    end
    
    alpha=0;
    
    for yy=1:param.maxiter
        alpha0=alpha;
        for i=1:dim
            v=partialprodscal(A,w,i);
            w.F{i}=reduce(v,i);
            w=normalizefuns(w);
            alpha=w.alpha;
            w.alpha=1;
        end
        erriter(yy)=abs(alpha-alpha0)/(alpha+alpha0);
        if param.display
            fprintf('  iteration #%d - stagnation = %d\n',yy,erriter(yy))
        end
        if erriter(yy)<param.itercrit
            break
        end
    end
    A=A-alpha*w;
    u=u+alpha*w;
    
    switch param.errorindicator
        case {'residual','reference'}
            errorder(zz)=norm(A)/normref;
        otherwise
            errorder(zz)= sqrt(min(u.alpha)^2/sum(u.alpha.^2));
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