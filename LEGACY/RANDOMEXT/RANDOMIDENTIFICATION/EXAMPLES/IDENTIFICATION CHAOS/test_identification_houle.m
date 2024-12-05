load('samples_H_T')
Xs = [Hs;Ts];

%%
I = randperm(size(Xs,2));
Xs = Xs(:,I(1:500));


%%
load houleidentifiee_chaos
%%
figure(100)
subplot(1,2,1)
multipdfsampleplot(Xs,'npts',30)
ax0 = axis;
cax0=caxis;
view(2)
xlabel('H')
ylabel('T')
title('Samples : Joint PDF H-T')
subplot(1,2,2)
multipdfplot(Xpc,'npts',30,'nbs',1e5,'axis',ax0)
axis(ax0)
caxis(cax0)
view(2)
xlabel('H')
ylabel('T')
title('Chaos decomposition : Joint PDF H-T')
%%
figure(101)
clf
pdfsampleplot(Xs(1,:),'npts',30,'ksdensity')
ax0 = axis;
hold on
pdfplot(Xpc(1),'k','nbs',1e4,'npts',30,'ksdensity','axis',ax0)
title('Marginal PDF H')
legend('samples','chaos decomposition')
figure(102)
clf
pdfsampleplot(Xs(2,:),'npts',30,'ksdensity')
ax0 = axis;
hold on
pdfplot(Xpc(2),'k','nbs',1e4,'npts',30,'ksdensity','axis',ax0)
title('Marginal PDF T')
legend('samples','chaos decomposition')

%%
pc = POLYCHAOS(2,3);

figure(12)
clf
subplot(1,2,1)
multipdfsampleplot(Xs,'npts',30)
view(2)
ax0 = axis;
cax0=caxis;
pause(2)

L=inf;
scanL=[];
scanXpc = cell(1,0);
%%
for i=1:50
    [Xpc0,L0] = PCidentification(pc,Xs,'nbestim',3e3,'rs',50,'fmin',20,'parametrized');
    if L0<L
        Xpc = Xpc0;
        L=L0;
        scanL=[scanL;L];
        scanXpc=[scanXpc , {Xpc}];
        scanL
        figure(12)
        subplot(1,2,2)
        multipdfplot(Xpc,'npts',30,'nbs',1e5,'axis',ax0)
        axis(ax0)
        caxis(cax0)
        pause(1)
    end
end



%% Courant

load -ascii currentgabonsurf

Xs = currentgabonsurf;
%%
figure(10)
clf
pdfsampleplot(Xs,'npts',60)
hold on
[P,x]=separatesamples(Xs,[32,50]);
[min(Xs),10,30,40,50,60,max(Xs)];
Xpc = eicdfpcproject(Xs,4,POLYFE(x),'nbgauss',30);

hold on
pdfplot(Xpc,'r','nbs',1e5,'npts',200)





%%
figure(10)
clf
pdfsampleplot(Xs,'npts',60)
hold on

leg = {'echantillons'};
scanp = [3,10,30,50];
for i=1:length(scanp)
    XpcV1D = eicdfpcproject(Xs,scanp(i),[],'nbgauss',60)
    figure(10)
    pdfplot(XpcV1D,getcourbestyles(i+1,'nomarker'),'nbs',1e5,'npts',100,'ksdensity')
    leg = [leg,{['chaos 1D ordre ' num2str(scanp(i))]}];
    pause(0.1);
end
legend(leg{:})
%%
figure(11)
clf
rep=find(Xs<=1|Xs>=69);Xstraf=Xs;
Xstraf(rep)=[];
pdfsampleplot(Xstraf,'npts',60)
hold on
leg = {'echantillons'};
scanp = [1,2,3,4,6,13];
scanp = 5;
for i=1:length(scanp)
    XpcV = eicdfpcproject_multimodes(Xstraf,scanp(i),POLYHERMITE(),[32,50],'nbgauss',30)
    figure(11)
%pdfplot(Xpc,'r','nbs',3e5,'npts',400)
    pdfplot(XpcV,getcourbestyles(i+1,'nomarker'),'nbs',1e5,'npts',100)
    leg = [leg,{['chaos 2D ordre ' num2str(scanp(i))]}];
end
legend(leg{:})
%%
XpcV = eicdfpcproject_multimodes(Xs,4,POLYHERMITE(),[10,30,40,50,60],'nbgauss',30)
figure(10)
pdfplot(XpcV,'r','nbs',2e5,'npts',200)
%%
[P,x,Xblock,mublock] = separatesamples(Xs,[32,50]);
pc = POLYCHAOS(RANDPOLYS(POLYFE(x),POLYHERMITE()),[0,3],'typebase',2);
[Xpc0,L0] = PCidentification(pc,Xs,'mublock',mublock,...
    'nbestim',3e3,'rs',10,'fmin',12,'parametrized','nocenter');
figure(10)
pdfplot(Xpc0,'g','nbs',3e5,'npts',400)

%%
XpcV = setnumber(XpcV,[3,4]);
PC = union(getPC(XpcHT),getPC(XpcV));
XpcHTV = [project(XpcHT,PC);project(XpcV,PC)];
PCu = POLYCHAOS(RANDPOLYS(PC),[6,6,0,6],'typebase',2)
PCu = restrictorder()
%%
figure(1)
clf
pdfsampleplot(currentgabonsurf,'npts',30)
ax0 = axis;
hold on
pdfplot(XpcHTV(3),'nbs',1e5,'r','npts',30,'ksdensity','axis',ax0)

figure(2)
clf
pdfsampleplot(Hs,'npts',30)
ax0 = axis;
hold on
pdfplot(XpcHT(1),'nbs',1e5,'r','npts',30,'ksdensity','axis',ax0)

figure(3)
clf
pdfsampleplot(Ts,'npts',30)
ax0 = axis;
hold on
pdfplot(XpcHTV(2),'nbs',1e5,'r','npts',30,'ksdensity','axis',ax0)



%% identification houle

Hpc = eicdfpcproject(Hs,15,POLYHERMITE(),'nbgauss',50)
figure(30)
clf
pdfsampleplot(Hs,'npts',40);
hold on
pdfplot(Hpc,'r','nbs',1e5,'npts',40)

%%

Xpc = eicdfpcproject2D(Xs,10,[],[],'nbgauss',30,'nbblock',15)
XpcHT = Xpc;
%%

Xpc = eicdfpcproject2D_multimodes([Ts;Hs],8,[],[],10.5,'nbgauss',30,'nbblock',15)
