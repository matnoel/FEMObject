function U=GSDthermiquesolvedetM(S,tol,ka,kn,kl,Ui,Urond,Utilde,Ubar,Uchap,varargin);

display_ = ischarin('display',varargin);

if display_
    fprintf('\n RESOLUTION PROB DETERMINISTE\n')
end

thermique_nonlin_forms

A = ka*a{S}(:,:);
if nbsec==1
    b =  kl*l{S}(:) ;
else
    b =  kl{1}*l{S}(:) + kl{2}*l2{S}(:) ;       
end

b = b - a{S}(:,Urond);
for i=1:length(Ui)
    A = A+2*g{S}(:,Ubar{i},Ui{i},:)+g{S}(:,:,Ubar{i},Ui{i});
    b = b - g{S}(:,Uchap{i,i},Ui{i},Ui{i});
    for j=i+1:length(Ui)
        b = b - 2*g{S}(:,Uchap{i,j},Ui{i},Ui{j});    
    end
end

U = zeros(getnbddlfree(S),1);

for k=1:100
    res = b - A*U - g{S}(:,Utilde,U,U) - 2*g{S}(:,U,Utilde,U)-kn*g{S}(:,U,U,U);
    nres = norm(res)/norm(b);

    if display_
        fprintf('iteration %d , residual = %.2d\n',k,nres);
    end
    if nres<tol
        break
    end

    AT = A+kn*(g{S}(:,:,U,U)+2*g{S}(:,U,U,:))+ ...
        2 * (g{S}(:,Utilde,:,U) + g{S}(:,U,Utilde,:)+ g{S}(:,:,U,Utilde));

%AT = A+kn*dn{S}(:,:,U,U);
%AT = A+kn*dn{S}(:,:,U,U)+ ...
%    2 * (g{S}(:,Utilde,:,U) + g{S}(:,U,Utilde,:)+ g{S}(:,:,U,Utilde));
%AT = A + dnapprox{S}(:,:,U,U) ;
%AT=A;
    dU = AT\res;
    U = U+dU;

end
if ~display_
    fprintf('(FM(l) %d iterations)',k);
end

if nres>=tol
    warning(['non convergence after #' num2str(k) ' iterations']);
end

