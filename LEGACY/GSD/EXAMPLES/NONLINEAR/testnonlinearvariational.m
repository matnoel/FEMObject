%% Definition probleme
r=15;
P = POINT([ 0,0 ; 1,0 ; 0,2 ; 1,2 ; 1,1 ; 2,1 ; 2,2 ; 0,1 ]);
S1 = mesh(DOMAIN(2,P(1),P(5)),r,r);
S2 = mesh(DOMAIN(2,P(8),P(4)),r,r);
S3 = mesh(DOMAIN(2,P(5),P(7)),r,r);
S = union(S1,S2);
S = union(S,S3);
S = concatgroupelem(S);

S = createddlnode(S,DDL('u'),DDL('g'));
S = addcl(S,create_boundary(S),'u',0);

RV = RANDVARS();
RV{1}=RVUNIFORM(0.7,1.3);
RV{2}=RVUNIFORM(0.7,1.3);
RV{3}=RVNORMAL(1,0.2);

X = PCMODEL(RV,'order',3,'pcg');


%%
a = BILINFORM(1,1)
l = LINFORM(0,1)

u = a{S}(:,:)\l{S}(:);
figure(1)
clf
plot(FENODEFIELD(u),S);

%%
Au = a{S}(:,u);
norm(Au-B)

%%
