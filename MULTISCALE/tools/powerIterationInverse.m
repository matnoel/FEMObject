function varargout = powerIterationInverse(A,v0,tol)
% function d = powerIterationInverse(A,v0,tol)
% function [v,d] = powerIterationInverse(A,v0,tol)
% function [v,d,flag] = powerIterationInverse(A,v0,tol)
% Power iteration inverse method for computing the inverse d of A's lowest magnitude
%   eigenvalue and a corresponding eigenvector v

d = 0;
u = v0;
v = u/norm(u);
residual = 1;
i = 0;
while residual > tol
    t = tic;
    i = i+1;
    d_old = d;
    u_old = u;
    if isfloat(A)
        u = A\v;
    else
        A = fcnchk(A);
        error('not implemented')
    end
    d = v'*u;
    v = u/norm(u);
    residual = abs(d-d_old)/abs(d);
    % residual = abs(norm(u)-norm(u_old))/norm(u);
    
    fprintf('\nIteration #%2.d : residual = %.3e w.r.t. eigenvalue',i,residual);
    fprintf('\n                eigenvalue = %f',d);
    fprintf('\n                elapsed time = %f s\n',toc(t));
end
fprintf('\n');
if nargout <= 1
    varargout{1} = d;
else
    varargout{1} = u(:);
    varargout{2} = d;
    if nargout == 3
        varargout{3} = double(residual<=tol);
    end
end

end
