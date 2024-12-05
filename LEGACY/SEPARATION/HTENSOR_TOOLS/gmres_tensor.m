function [x, norm_r] = gmres_tensor(apply_mat, apply_precond, b, opts, x0)
% Truncated GMRES method for htensor.
%
%   [X, NORM_R] = GMRES_TENSOR(APPLY_MAT, APPLY_PRECOND, B, OPTS) applies the
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
%   - OPTS.RESTART     maximal restart
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
norm_b = norm(b);

v = cell(opts.restart,1);
w = v;

for i=1:opts.maxit
    r = b-apply_mat(x,opts);
    norm_r(i) = norm(r)/norm_b;
    fprintf('GMRES(m) : Iter #%d - Relative Residual %d\n',...
        i,norm_r(i));
    if(opts.plot_conv)
        semilogy(norm_r, 'b');
        hold on;
        drawnow;
    end
    if( norm_r(i) < opts.tol || ...
            (isfield(opts, 'rel_tol') && (i>1) && ...
            norm_r(i-1) - norm_r(i) < opts.rel_tol*norm_r(i)) )
        return
    end
    v{1} = truncate(apply_precond(r,opts),opts);
    v{1} = (1/norm(v{1})) * v{1};
    for j=1:opts.restart
        w{j} = apply_precond(apply_mat(v{j},opts),opts);
        w{j} = (1/norm(w{j}))*w{j};
        VV = buildVV(v,j);
        Vw = buildVw(v,w{j},j);
        a = VV\Vw;
        v{j+1} = w{j};
        for k = 1:j
            v{j+1} = v{j+1}-a(k)*v{k};
        end
        v{j+1} = truncate(v{j+1},opts);
    end
    WAV = buildWAV(w,v,apply_mat,opts);
    Wr = buildWr(w,r,opts);
    y = WAV\Wr;
    for j=1:opts.restart
        x = x + y(j)*v{j};
    end
    x = truncate(x,opts);
end

end
%%

function VV = buildVV(v,k)
VV = zeros(k);
for i=1:k
    for j=1:k
        VV(i,j) = innerprod(v{i},v{j});
    end
end
end

function Vw = buildVw(v,w,k)
Vw = zeros(k,1);
for i=1:k
    Vw(i) = innerprod(v{i},w);
end
end

function WAV = buildWAV(w,v,apply_mat,opts)
WAV = zeros(opts.restart);
for i=1:opts.restart
    for j=1:opts.restart
        WAV(i,j) = innerprod(w{i},apply_mat(v{j},opts));
    end
end
end

function Wr = buildWr(w,r,opts)
Wr = zeros(opts.restart,1);
for i=1:opts.restart
    Wr(i) = innerprod(w{i},r);
end
end