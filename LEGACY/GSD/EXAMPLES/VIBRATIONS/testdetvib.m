%% modele elements finis (S, f, K1, M1, ddloutput)
choixmesh = 1;

switch choixmesh
    case 1
        example_2Dsimple_gmsh
    case 2
        example_piece_gmsh
    case 3
        example_plaque_gmsh
end
E = 1;
RHO = 1;

M = RHO*M1;
K = E*K1;

% DEFINITION DU MODELE 

alphaw = 0.01;
betaw = 0.01;
wband = [0,1];
cout = 2;

%%
scanw = linspace(wband(1),wband(2),20);

[urw,uiw,umodw,uphiw] = calc_frf(scanw,K,alphaw*M+betaw*K,M,f,ddloutput);

%

figure(100)
clf
subplot(1,2,1)
semilogy(scanw,umodw(cout,:),'r-')
hold on
subplot(1,2,2)
plot(scanw,uphiw(cout,:),'r-')
hold on


%%
p=4;
[ur,ui,upc]=calc_frf_spectral(linspace(wband(1),wband(2),10),p,K,M,f,alphaw,betaw);

%% RESOLUTION GSD
GSD = GSDSOLVER('tol',1e-2,'nbfoncmax',20,...
    'display',true,'update',true,'errorindicator','residual',...
    'inittype','random','direct',true,'reuse',true);
GSDA = setparam(GSD,'type','arnoldi','restart',3,'orthocrit',1e-11,'nbfoncmaxsimul',20);
GSDP = setparam(GSD,'type','power','pfixmax',3,'pfixtol',1e-3);
GSD=GSDP;

p=4;
%[urad] = solve(GSDA,Aw,bw);
[ur,ui,urad]=calc_frf_spectral(linspace(wband(1),wband(2),10),p,K,M,f,alphaw,betaw,@(A,b) solve(GSDP,A,b));


%%
p=3;
crit = 1e-2;
xmesh0 = [wband(1),wband(2)/3;wband(2)/3,2*wband(2)/3;2*wband(2)/3,wband(2)];

figure(12)
clf
uh=[];
for r=1:10
    xmesh0
    xmesh=xmesh0;
    n=size(xmesh,1);
    GSD = setparam(GSD,'nbfoncmax',10);
    [ur,ui,uh]=calc_frf_spectral(xmesh,p,K,M,f,alphaw,betaw,@(A,b) solve(GSD,A,b,uh));
    PC = getPC(uh);
    GSD = setparam(GSD,'nbfoncmax',0);
    [ur,ui,uhfin]=calc_frf_spectral(xmesh,p+1,K,M,f,alphaw,betaw,@(A,b) solve(GSD,A,b,uh));
    PCfin = getPC(uhfin);
    uh = project(uh,PCfin);
    ehd = double(expand(uh(ddloutput(cout))-uhfin(ddloutput(cout))));

    repelem = mod(0:length(PCfin)-1,n)+1;
    errorelem=zeros(1,n);
    splitelem=zeros(1,n);
    xmesh0 = zeros(0,2);
    uhfind = double(expand(uhfin(ddloutput(cout))));

    for k=1:n 
        errorelem(k) = sqrt(sum(sum(ehd(:,repelem==k).^2))./sum(sum(uhfind(:,repelem==k).^2)));
        splitelem(k) =  errorelem(k)>crit;
        if splitelem(k)
            bar(sum(xmesh(k,:),2)/2,errorelem(k),xmesh(k,2)-xmesh(k,1),'r');
            xmid = sum(xmesh(k,:))/2;
            xmesh0 = [xmesh0;[xmesh(k,1),xmid];[xmid,xmesh(k,2)]];    
        else
            bar(sum(xmesh(k,:),2)/2,errorelem(k),xmesh(k,2)-xmesh(k,1),'w');
            xmesh0 = [xmesh0;xmesh(k,:)];   
        end 

        hold on
    end

    plot(unique(xmesh),zeros(size(unique(xmesh))),'.')
    pause(1)
    if ~any(splitelem)
        break
    end

    urad=uh;
    hold off
end


%%
scanwrad = linspace(wband(1),wband(2),100);
uradr = urad(1:size(K,1));
uradrout = uradr(ddloutput);
uradi = urad(size(K,1)+1:end);
uradiout = uradi(ddloutput);
clear j
scanwrad = linspace(wband(1),wband(2),300);
uradroutw = double(randomeval(uradrout,scanwrad',RANDVARS(RVUNIFORM(wband(1),wband(end)))));
uradioutw = double(randomeval(uradiout,scanwrad',RANDVARS(RVUNIFORM(wband(1),wband(end)))));
uradmodw = abs(uradroutw+j*uradioutw);
uradphiw = angle(uradroutw+j*uradioutw);


%

figure(100)
subplot(1,2,1)
semilogy(scanwrad,uradmodw(cout,:),'b-')
hold on
subplot(1,2,2)
plot(scanwrad,uradphiw(cout,:),'b-')
hold on
