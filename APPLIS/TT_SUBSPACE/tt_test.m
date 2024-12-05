clear all
%% Basic setup
d = 10;
n = 100;
pb = tt_poisson(d,n+1:n+d+1);
A = pb.Att; % A.core is directly converted from canonical to TT
            % format, we know that we can do better
b = pb.btt;

%%
opts.maxiter = 10;
opts.errorindicator = 'residual';
u0 = greedyapproxsol(A,b,[],opts);
u0 = LRTENSOR(TT_CORE(u0.core),u0.space);
%%
u = orth(u0);

stag = @(x0,x) norm(orth(x-x0))/norm(orth(x0));
res = @(x) norm(orth(A*x-b))/norm(orth(b));

u_als = tt_als(A,b,u);
u_dmrg = tt_dmrg(A,b,u);

stag_als = stag(u,u_als)
stag_dmrg = stag(u,u_dmrg)
res_als = res(u_als)
res_dmrg = res(u_dmrg)

%% First version of the algorithm
opts.optim_core = 'als';
opts.max_iter = 10;
u1 = tt_subspace(A,b,opts);

%% 
opts.optim_core = 'dmrg';
opts.max_iter = 10;
u2 = tt_subspace(A,b,opts);



