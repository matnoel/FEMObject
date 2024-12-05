

adrot_example

%%
nXi=100;
Xi = linspace(0,1,nXi)';



%% Operateur A

AT = TSPACE_OPERATORS( [{A(:)} ;{param_operator(PhiA(Xi))} ] );
AC = CANONICAL_CORE( ones(1,AT.dim(1)),2 );
A_LRT = LRTENSOR(AC,AT)

%% Right hand side

BT = TSPACE_VECTORS( [{b} ; {ones(nXi,1)} ] );
BC = CANONICAL_CORE( ones(1,BT.dim(1)),2 );
B_LRT = LRTENSOR(BC,BT)

%% Inverse interpolation

XiK = floor( ([0.05 0.2 0.8]')*nXi ) ;
Ym = eval_param_operator(A_LRT,XiK);
Ym = cellfun(@(y)implicit_inverse_lu(y),Ym,'UniformOutput',0);

%% Frob projection
[P,residual] = inverse_interp_frob(A_LRT,Ym)

%% Semi Frob. projection

H=my_hadamard(A_LRT.sz(1),20);
[P_H,residual_H] = inverse_interp_frob(A_LRT,Ym,H)

%% Compare the interpolation function 
tp1=cellfun( @spdiags,P.space.u{2},'UniformOutput',0);
tp1=[tp1{:}];
tp2=cellfun( @spdiags,P_H.space.u{2},'UniformOutput',0);
tp2=[tp2{:}];
plot(tp1)
hold on
plot(tp2,'--')
hold off


%% Greedy frob ?

m=3;

[P,residual]=greedy_inverse_interp_frob(A_LRT,Ym,m,H);







%%





