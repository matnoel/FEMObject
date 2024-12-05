function [kappa,beta,alpha] = condest_PA(P,A,k_eval)

%% LAISSE TEBON .... !!!!



% Give an estimation of the condition number kappa(PA) at the point k
% Compute the largest (resp.smallest) singular value using (inverse) power
% method.

A_eval = eval_sparse(A,k_eval);

PA = @(v) cell2mat(cellfun( @(L,U) U\(L\(A_eval*v)) ,P.L,P.U,'UniformOutput',0)') * P.Phi(:,k_eval) ;
PAt= @(v) A_eval'* ( cell2mat(cellfun( @(L,U) ( (v'/U) /L)' ,P.L,P.U,'UniformOutput',0)') * P.Phi(:,k_eval) );

P_fun = @(v) cell2mat(cellfun( @(L,U) U\(L\v) ,P.L,P.U,'UniformOutput',0)') * P.Phi(:,k_eval) ;
Pt_fun= @(v) ( cell2mat(cellfun( @(L,U) ( (v'/U) /L)' ,P.L,P.U,'UniformOutput',0)') * P.Phi(:,k_eval) );

% PA = @(v) eval(P*(A*v) ,k) ;
% PAt= @(v) A'* ( cell2mat(cellfun( @(L,U) ( (v'/U) /L)' ,P.L,P.U,'UniformOutput',0)') * P.Phi(:,k_eval) );





K=500;
tol = 1e-4;

v=randn(size(A_eval,1),1);
v=v/norm(v);
k=0;
stag=100;
while (k<K) && (stag>tol)
    
    k=k+1;
    
    vtmp = P_fun(A_eval*v);
    vtmp = A_eval'* Pt_fun(vtmp);
    
    beta = norm(vtmp);
    vtmp = vtmp/beta;
    stag = norm(v-vtmp);
    
    v=vtmp;
end

beta=sqrt(beta);




%%

tol_cgs=1e-6;



tmp=cellfun( @(L,U) U\(L\full(A_eval)) ,P.L,P.U,'UniformOutput',0);
PA_full=tmp{1}*P.Phi(1,k_eval) ;
for k=2:P.r
    PA_full=PA_full+tmp{k}*P.Phi(k,k_eval) ;
end


v=randn(size(A_eval,1),1);
v=v/norm(v);
k=0;
stag=100;


[sol,flag,relres,iter,resvec] = gmres( PA_full  ,v ,15,tol_cgs )



% Pm  = @(v) funalacon( @(x) cgs( P_fun  ,x ,tol_cgs ),v );
% PAmt = @(v) funalacon( @(x) cgs( Pt_fun ,x ,tol_cgs ),v );

PAm  = @(v) funalacon( @(x) gmres( PA  ,x ,15,tol_cgs ),v );
PAmt = @(v) funalacon( @(x) gmres( PAt ,x ,15,tol_cgs ),v );



% PAm  = @(v) funalacon( @(x) cgs( PA   ,x ,tol_cgs ),v );
% PAmt = @(v) funalacon( @(x) cgs( PAt  ,x ,tol_cgs ),v );


%%% avec l'affichage cgs
% PAm  = @(x) cgs( PA  ,x ,tol_cgs );
% PAmt = @(x) cgs( PAt ,x ,tol_cgs );



% PAm  = @(v) funalacon( @(x) cgs( @(w) P_fun(A_eval*w)   ,x ,tol_cgs ),v );
% PAmt = @(v) funalacon( @(x) cgs( @(w) A_eval'*Pt_fun(w) ,x ,tol_cgs ),v );
% 
% 
% PAm  = @(x) cgs( @(v) P_fun(A_eval*v)   ,x ,tol_cgs );
% PAmt = @(x) cgs( @(v) A_eval'*Pt_fun(v) ,x ,tol_cgs );

v=randn(size(A_eval,1),1);
v=v/norm(v);
k=0;
stag=100;


while (k<K) && (stag>tol)
    
    k=k+1;
    
    vtmp = PAmt(v);
    vtmp = PAm(vtmp);
    
    alpha=norm(vtmp);
    vtmp=vtmp/alpha;
    stag= norm(v-vtmp)
    
    v=vtmp;
end

alpha=1/sqrt(alpha);


%%

kappa=beta/alpha;


end

function [a,b]=funalacon(v,x)
[a,b]=v(x);
end