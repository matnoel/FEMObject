%% Construction d'une base POD :

KMAX=100;
r = 10;

U_POD =[];
for k=1:KMAX
    ind=floor( rand()*nXi )+1;
    U_POD = [U_POD , eval_sparse(A_Xi,ind) \ b ];
end

[U_POD,~,~]=svd(U_POD,'econ');
U_POD=U_POD(:,1:r);


%% projection de galerkin

UPAU = (Precond'*U_POD)' * (A_Xi*U_POD);

UPAU = compress(UPAU);
UPb  = (Precond'*U_POD)'*b;
UAU  = ((A_Xi*U_POD)'*U_POD)';

%%
k=80;

u_g = U_POD* ( eval(UAU,k)\(U_POD'*b ) );
u_r = U_POD*(eval(UPAU,k)\eval(UPb,k));
u_e = eval_sparse(A_Xi,k)\b;
u_p = U_POD*(U_POD'*u_e);

disp('meilleure approx')
disp( norm(u_p-u_e)/norm(u_e) )
disp('Galerkin')
disp(norm(u_g-u_e)/norm(u_e))
disp('Galerkin precond')
disp( norm(u_r-u_e)/norm(u_e) )


%% Reduced basis

rank=10;
[U,maxNPr,xie,ResNorm,error_P]=greedy_preconditioned_RB(A_Xi,b,Precond,rank,Xi)

cellfun(@max,error_P)

%%

[U,xie,error_E]=greedy_exact(A_Xi,b,rank,Xi);
error_E
%%


I = affine_matrix( {speye(A_Xi.s)}, ones(1,A_Xi.n) );
[U,maxNPr,xie,ResNorm,error_I]=greedy_preconditioned_RB(A_Xi,b,I,rank,Xi);


cellfun(@max,ResNorm)







