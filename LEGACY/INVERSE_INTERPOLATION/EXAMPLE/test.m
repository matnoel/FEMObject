

nXi=100;

% Xi = lhsdesign(nXi,dim_p);
Xi = linspace(0,1,nXi)';

A_Xi = affine_matrix( A(:),PhiA(Xi));


%%

Xi_interp = [15 30 45 60 85 100];

A_Xi_interp={};
for k=1:length(Xi_interp)
    A_Xi_interp{k} = eval_sparse(A_Xi,Xi_interp(k));
end

Precond = implicit_affine_matrix(A_Xi_interp,[]);

[Precond,M,S,out] = frobenius_projection_exacte(Precond,A_Xi,1,0);

%% Frob semi norm

dimH=10;

% H = my_hadamard(A_Xi.s(1), dimH );
H = 2*(rand(A_Xi.s(1),dimH)<0.5)-1;


[Precond,M,S,out] = frobenius_projection_V(Precond,A_Xi,H,1,0);

%%
plot(Precond.Phi')
hold on
for k=1:length(Xi_interp)
    plot( [Xi_interp(k) Xi_interp(k)],[0 1] ,'--+k')
end
hold off

%%

kappa=compute_condition_number(A_Xi,Precond)
figure(2)
semilogy(kappa)

%%

















