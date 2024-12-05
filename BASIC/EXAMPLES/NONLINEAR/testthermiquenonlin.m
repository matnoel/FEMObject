r=10;
P = POINT([ 0,0 ; 1,0 ; 0,2 ; 1,2 ; 1,1 ; 2,1 ; 2,2 ]);
S1 = mesh(DOMAIN(2,P(1),P(4)),r,r*2);
S2 = mesh(DOMAIN(2,P(5),P(7)),r,r);
S = union(S1,S2);

mat=FOUR_ISOT('k',1,'k2',1);
S = setmaterial(S,mat);
S=final(S);

S=addcl(S,create_boundary(S),'T',0);
scanf=[3,5,20];
figure(1)
clf

N = NEWTONSOLVER('type','modified','tol',1e-6,'tolreact',1e-1);

for kk=1:length(scanf)

    f=bodyload(S,[],'QN',scanf(kk));

    K=calc_rigi(S);

    qlin=K\f;

    q = solve(N,f,@(u) calc_fint(S,u),@(u) calc_rigitang(S,u));

    q=unfreevector(S,q);
    qlin=unfreevector(S,qlin);

    figure(1)
    subplot(2,length(scanf),kk)
    plot(FENODEFIELD(qlin),S,'surface')
    c0 = caxis;
    a0=axis;
    subplot(2,length(scanf),kk+length(scanf))
    plot(FENODEFIELD(q),S,'surface')
    caxis(c0);
    axis(a0)

end


