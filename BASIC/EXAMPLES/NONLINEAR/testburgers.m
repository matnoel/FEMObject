figure(1)

scanr = [26];
for r=scanr
    x0=-10;
    x1=10;
    meshtype=2;
    P = POINT([x0;x1]);
    switch meshtype 
    case 1

        S = mesh(DOMAIN(1,P(1),P(2)),r);

    case 2
        % funmesh = inline('((x-1/2).^3+a*(x-1/2)+a/2+1/8)/2/(1/8+a/2)','x','a');
        funmesh = @(x,a) ((x-1/2).^3+a*(x-1/2)+a/2+1/8)/2/(1/8+a/2);

        xx = linspace(0,1,r+1);
        xx = x0 + funmesh(xx,1/8)*(x1-x0);

        N = NODE(POINT(xx(:)));
        S = MODEL('UNID');
        S = addnode(S,N);
        E = SEG2(N,1:numel(N)-1,[[1:numel(N)-1]',[2:numel(N)]']);
        S = addelem(S,E);
    end
    figure(16)
    clf
    plot(S,'node')


    scank=[0.1:0.1:0.5];
    for i=1:length(scank)
        k=scank(i)
        mat=BURGERS('k',k,'formulation',1,'alpha',0);
        S = setmaterial(S,mat);
        S=final(S);
        mu0 = 1;
        S=addcl(S,P(1),'u',1/2*(1+tanh(x0/4/mu0)));
        S=addcl(S,P(2),'u',1/2*(1+tanh(x1/4/mu0)));
        u0 = calc_init_dirichlet(S);


        N = NEWTONSOLVER('type','full','tol',1e-6,'tolreact',1e-1,'maxiter',20);

        f=freevector(S,zeros(S.nbddl,1));
        q = solve(N,0,@(u) calc_fint(S,unfreevector(S,u)),@(u) calc_rigitang(S,unfreevector(S,u)),u0,[],[],S);
        u0 = q ; 
        figure(1)
        plot(getcoord(S.node),q,'g');
        hold on
        pause(.2)


    end

end
