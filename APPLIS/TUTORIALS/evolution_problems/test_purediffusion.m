
S=gmsh(DOMAIN(2),0.03);
xnode = getcoord(getnode(S));
figure(1)
clf
plot(S)

%%
casmu = 0

switch casmu
case 0
    k = 1; 
    mutext = '';    
case 1
    k = 1e-2; 
    mutext = 'mu1e-2_';
case 2
    k = 1e-4; 
    mutext = 'mu1e-4_';
case 3
    k = 1e2; 
    mutext = 'mu1e2_';                
case 4
    k = 1e4; 
    mutext = 'mu1e4_';  
case 5
    k = 1e-6; 
    mutext = 'mu1e-6_';
case 6
    k = 1e6; 
    mutext = 'mu1e6_';
end

manu=0;

mat = FOUR_ISOT('k',1,'c',1);
S = setmaterial(S,mat);
S=final(S);
S=addcl(S,[],'T');
Mx = calc_freematrix(S,@mass);
A = k*calc_freematrix(S,@diff);
N = DGTIMESOLVER(TIMEMODEL(0,1,100),0);


if manu==1
    x = xnode(:,1);y=xnode(:,2);
    t = gettapprox(N);
    factors = [1,1/2,1/2,1/3,1/5];
%factors = [1,1,1,1,1];
    uref = factors(1)*(sin(pi*x).*sin(pi*y)).*N(@(t) t)+...
        factors(2)*(sin(pi*x).*sin(2*pi*y)).*N(@(t) (1-t))+...
        factors(3)*(sin(2*pi*x).*sin(pi*y)).*N(@(t) sin(pi*t))+...
        factors(4)*(sin(2*pi*x).*sin(2*pi*y)).*N(@(t) sin(2*pi*t))+...
        factors(5)*(sin(4*pi*x).*sin(4*pi*y)).*N(@(t) sin(4*pi*t));
    uref = freevector(uref,S);
    b = Mx*diff(uref)+A*uref;
    fichname = 'example_diff_manu_' ;


else
    b = bodyload(S,[],'QN',1);
    bbis = bodyload(S,[],'QN',@(x) x(1));
    b = b*one(N) + 2*bbis*rampe(N,0,gett1(N));
    uref = dsolve(N,b,Mx,A);
    fichname = ['example_diff_' mutext] ;
end
%

figure(8)
clf
evol(uref,S,'surface','rescale','z')
%evol(ut,S)


save(fichname,'A','b','Mx','uref','N','S')


Mxmetric = A;
Mtmetric = getMmatrix(N);
