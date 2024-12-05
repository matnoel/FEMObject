function [kappa,beta,alpha]=my_cond_precond(M,PP)


K=500;
tol = 1e-2;

v=randn(size(M,1),1);
v=v/norm(v);
k=0;
stag=100;
while (k<K) && (stag>tol)
    
    k=k+1;
    
    vtmp=PP*(M*v);
    vtmp=(M')*((PP')*vtmp);
    
    beta=norm(vtmp);
    vtmp=vtmp/beta;
    stag= norm(v-vtmp);
    
    v=vtmp;
end

beta=sqrt(beta);

%%

[LM,UM]=lu(M);
[LP,UP]=lu(PP);

v=randn(size(M,1),1);
v=v/norm(v);
k=0;
stag=100;


while (k<K) && (stag>tol)
    
    k=k+1;
    
    vtmp= UP\(LP\v);
    vtmp= UM\(LM\vtmp);
    
    vtmp= (LM')\((UM')\vtmp);
    vtmp= (LP')\((UP')\vtmp);
    
    alpha=norm(vtmp);
    vtmp=vtmp/alpha;
    stag= norm(v-vtmp);
    
    v=vtmp;
end
alpha=1/sqrt(alpha);

%%

kappa=beta/alpha;








