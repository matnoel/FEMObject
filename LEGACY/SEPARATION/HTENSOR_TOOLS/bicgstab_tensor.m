function [x, norm_r] = bicgstab_tensor(apply_mat, apply_precond, b, opts, x0)
%Truncated BICGSTAB method for htensor.
%
%   [X, NORM_R] = BICGSTAB_TENSOR(APPLY_MAT, APPLY_PRECOND, B, OPTS) applies the
%   BICGSTAB method to solve the linear system
%      APPLY_MAT(X) = B,
%   using the preconditioner APPLY_PRECOND. The right-hand side B as well
%   as the solution are htensor objects. The iterates of BICGSTAB are truncated
%   to low hierarchical rank in each iteration, as specified by OPTS.
%
%   APPLY_MAT(X,OPTS) and APPLY_PRECOND(X,OPTS) are function handles that
%   accept an htensor input X and return a matrix-vector product as
%   htensor.
%
%   Required fields of OPTS:
%   - OPTS.MAXIT       maximal number of iterations of CG method
%   - OPTS.TOL         tolerance for residual norm in CG method
%   - OPTS.MAX_RANK    argument to TRUNCATE
%
%   Optional fields of OPTS:
%   - OPTS.ABS_EPS     argument to TRUNCATE
%   - OPTS.PLOT_CONV   if true, residual norm is plotted during iteration
%   - OPTS.REL_EPS     argument to TRUNCATE
%
%   See also APPLY_LIN_MAT, APPLY_INV_MAT.

if(~isfield(opts, 'plot_conv'))
  opts.plot_conv = true;
end

if(nargin == 5)
  x = x0;
else
    if all(rank(b)==1)
        x = 0*b;
    else
        x = truncate(0*b, struct('max_rank', 1));
    end
end

norm_r = zeros(1,opts.maxit);

r = b - apply_mat(x,opts);
rt = r;
rho0 = innerprod(rt,r);

for i = 1:opts.maxit
    rho = innerprod(rt,r);
    if rho == 0
        error('Method fails');
    end
    if i == 1
        p = r;
    else
        beta = (rho/rho0)*(alp/omega);
        rho0 = rho;
        p = r + beta*(p-omega * v);
    end
    p = truncate(p,opts);
    ph = apply_precond(p,opts);
    ph = truncate(ph,opts);
    v = apply_mat(ph,opts);
    v = truncate(v,opts);
    alp = rho / innerprod(rt,v);
    s = r - alp*v;
    s = truncate(s,opts);
    if norm(s) < opts.tol
        x = x + alp*p;
        return
    end
    sh = apply_precond(s,opts);
    sh = truncate(sh,opts);
    t = apply_mat(sh,opts);
    t = truncate(t,opts);
    omega = innerprod(t,s)/innerprod(t,t);
    x = x + alp*ph+omega*sh;
    x = truncate(x,opts);
    r = b - apply_mat(x,opts);
    r = truncate(r,opts);
    norm_r(i) = norm(r)/norm(b);
    fprintf('BICGSTAB : Iter #%d - Relative Residual %d\n',...
        i,norm_r(i));
    if(opts.plot_conv)
        semilogy(norm_r, 'b');
        hold on;
        drawnow;
    end
    if( norm_r(i) < opts.tol || ...
            (isfield(opts, 'rel_tol') && (i>1) && ...
            norm_r(i-1) - norm_r(i) < opts.rel_tol*norm_r(i)) )
        break;
    end
    if norm_r(i) > 1e3
        break;
    end
    if omega == 0
        error('Method fails')
    end
end

