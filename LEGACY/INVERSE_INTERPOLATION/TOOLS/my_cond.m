function [kappa,beta,alpha]=my_cond(M)


K=500;
tol = 3e-3;

v=randn(size(M,1),1);
v=v/norm(v);
k=0;
stag=100;
while (k<K) && (stag>tol)
    
    k=k+1;
    
    vtmp=M*v;
    vtmp=M'*vtmp;
    
    beta=norm(vtmp);
    vtmp=vtmp/beta;
    stag= norm(v-vtmp);
    
    v=vtmp;
end

beta=sqrt(beta);

%%

[L,U]=lu(M);


v=randn(size(M,1),1);
v=v/norm(v);
k=0;
stag=100;


while (k<K) && (stag>tol)
    
    k=k+1;
    
    vtmp= U\(L\v);
    vtmp= (L')\((U')\vtmp);
    
    alpha=norm(vtmp);
    vtmp=vtmp/alpha;
    stag= norm(v-vtmp);
    
    v=vtmp;
end
alpha=1/sqrt(alpha);

%%

kappa=beta/alpha;








