caseload = 1;

cx = @(x) (x(:,2)-1/2);
cy = @(x) (1/2-x(:,1));

S = gmsh(DOMAIN(2),0.015);


cx = @(x) (x(:,2)-1/2);
cy = @(x) (1/2-x(:,1));

S = gmsh(DOMAIN(2),0.015);


xnode = getcoord(getnode(S));
CX = cx(xnode);
CY = cy(xnode);
Cnorm = sqrt(sum(CX.^2+CY.^2,2));
C = [CX,CY];
Cnormalize = [CX./Cnorm,CY./Cnorm];

S=createddlnode(S,DDL('U'));
S=addcl(S,[],'U');

%%
caseload = 1;


amass = BILINFORM(0,0);
aeffectdiff = BILINFORM(1,1);
adisp1 = setfree(DIFFUSIONFORM(Cnorm,Cnormalize,Cnormalize),1);
adisp2 = setfree(DIFFUSIONFORM(Cnorm),1);
aadv = setfree(ADVECTIONFORM(C),1);

Amass = amass{S}(:,:);
Aeffectdiff = aeffectdiff{S}(:,:);
Adisp1 = adisp1{S}(:,:);
Adisp2 = adisp2{S}(:,:);
Aadv = aadv{S}(:,:);

funeval = @(u,x) transpose(reshape(eval_sol(S,unfreevector(S,full(double(expand(u)))),POINT(x),'T'),getnbtimedof(u),size(x,1)));
n=getnbddl(S);

figure(100)
clf
plot(S)

%
if caseload==1
    lambda = .00;
    de = .0001;
    alphal =0.1;
    alphat =.0001;
    ct = 3*pi;
else
    lambda = .0000;
    de = .00001;
    alphal =0.1;
    alphat =.0001;
    ct = 1;
end

astab = setfree(STABFORM(ct*C,1,de),1);
Astab =  astab{S}(:,:);


T = TIMEMODEL(0,1,100);
N = DGTIMESOLVER(T,1,'lu',true);
N = setevolparam(N,'plotstep',1);
N = setparam(N,'display',true);
N = setparam(N,'rescale','z');
t = gettapprox(N);


Mx = Amass;
A = lambda*Amass + ct*Aadv+ de*Aeffectdiff + ...
    ct*alphat*(Adisp2-Adisp1) + ...
    ct*alphal*Adisp1 + ...
    Astab;

switch caseload
case 0
    fsecform = LINFORM(0);
    alpha = 0.07;
    u0fun = @(x) exp(-((x(:,1)-2/3).^2+(x(:,2)-2/3).^2)/alpha^2);
    u0 = freevector(S,u0fun(xnode));figure(66);plot(u0,S,'surface')
    b = (Mx*u0)*getinitialequivalent(N);
    foutputfree = fsecform{S}(:);
    foutput = unfreevector(S,foutputfree)*one(N);
case 1
%fsec = fsecform{keepgroupelem(S,2)}(:);
    alpha = 0.07;
    f0fun = @(x) exp(-((x(:,1)-3/4).^2+(x(:,2)-1/2).^2)/alpha^2);
    fsecform = LINFORM(0,f0fun(xnode),0);
    fsec = fsecform{S}(:);
    funt = N(@(t) max((1-t/0.01),0));
    funt = N(@(t) max((1-t/0.01)>=0,0));
    funt=one(N);
    b = fsec*funt;
    foutputfree = fsecform{keepgroupelem(S,3)}(:);
    foutput = unfreevector(S,foutputfree)*one(N);
end



uref = dsolve(N,b,Mx,A);

figure(4);
clf
evol(uref,S,'rescale','z')


fichname = 'example_disp_' ;

save(fichname,'A','b','Mx','uref','N','S')


Mxmetric = Mx;%Aeffectdiff;
Mtmetric = getMmatrix(N);
