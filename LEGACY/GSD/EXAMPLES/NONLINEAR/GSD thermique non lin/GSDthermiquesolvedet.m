function U=GSDthermiquesolvedet(S,tol,ka,kn,kl,varargin)

display_ = ischarin('display',varargin);

if display_
    fprintf('\nRESOLUTION PROB DETERMINISTE\n')
end

thermique_nonlin_forms

%a = setfact(a,ka);
%n = setfact(n,kn);
%dn = setfact(dn,kn);
%l = setfact(l,kl);
%g = setfact(g,kn);

U0 = zeros(getnbddlfree(S),1);
U=U0;

A = ka*a{S}(:,:);
if nbsec==1
    b =  kl*l{S}(:) ;
else
    b =  kl{1}*l{S}(:) + kl{2}*l2{S}(:) ;    
end
inc = 0;
for k=1:100
    res = b - ka*a{S}(:,U) - kn*g{S}(:,U,U,U);

    nres = norm(res)/norm(b);
    if display_
        fprintf('iteration %d , residual = %.2d\n',k,nres);
    end
    if nres<tol
        break
    end

%AT = A+kn*dn{S}(:,:,U,U);
    AT = A+kn*g{S}(:,:,U,U)+2*kn*g{S}(:,U,:,U);
%AT = A + kn*dnapprox{S}(:,:,U,U);
%AT=A;
    dU = AT\res;
    U = U+dU;

%norm(U-U0)/norm(U)

    U0=U;
end

if display_
    fprintf('(F(l) %d iterations)',k);
end


if nres>=tol
    warning('NON CONVERGENCE');
end
