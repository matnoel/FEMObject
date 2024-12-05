
%% IDENTIFICATION 1 VARIABLE ALEATOIRE (OPTIM SUR HYPERSPHERE)
%----------------------------------------------------------------------------
%----------------------------------------------------------------------------
%----------------------------------------------------------------------------
%----------------------------------------------------------------------------

X = RVLOGNORMAL(1,0.3,0,'stat');
Xs = random(X,1,300);

pc = POLYCHAOS(1,2);

%% Vision non parametree
% l'algorithme fmin sous contrainte ne semble pas fonctionner (a voir)
% l'algo condor marche mieux (a voir pour la tolerance de la contraintes)
Xpc = PCidentification(pc,Xs,'plotfun',30,'rs',0,'parametrized');
Xpc = PCidentification(pc,Xs,'nbestim',3e3,'rs',2,'fmin',12,'illustrate','plotfun',10);

%% Vision parametree

Xpc = PCidentification(pc,Xs,'plotfun',30,'rs',0,'parametrized');
Xpc = PCidentification(pc,Xs,'nbestim',3e3,'rs',10,'fmin',10,'illustrate','parametrized');

%% Vision ou on n'impose pas la moyenne
Xpc = PCidentification(pc,Xs,'plotfun',15,'rs',0,'parametrized','nocenter');
Xpc = PCidentification(pc,Xs,'nbestim',3e3,'rs',10,'fmin',10,'illustrate','parametrized','nocenter');

%% IDENTIFICATION 2 VARIABLES ALEATOIRES 
%----------------------------------------------------------------------------
%----------------------------------------------------------------------------
%----------------------------------------------------------------------------
%----------------------------------------------------------------------------
X = RANDVARS(RVNORMAL(0,1),RVNORMAL(0,1));
X = get(PCMODEL(X,'order',1),'X');
%X = [X(1)+3*X(2),X(1)-2*X(2)];
Xs = double(full(random(X,1,500)));

pc = POLYCHAOS(2,1);

Xpc = PCidentification(pc,Xs,'plotfun',30,'rs',0,'parametrized');
Xpc = PCidentification(pc,Xs,'nbestim',3e3,'rs',10,'fmin',10,'illustrate','parametrized');


%% IDENTIFICATION 2 VARIABLES ALEATOIRES 
%----------------------------------------------------------------------------
%----------------------------------------------------------------------------
%----------------------------------------------------------------------------
%----------------------------------------------------------------------------
X = RANDVARS(RVLOGNORMAL(1,.5,0,'stat'),RVLOGNORMAL(1,.5,0,'stat'));
X = get(PCMODEL(X,'order',2),'X');
%X = [X(1)+3*X(2),X(1)-2*X(2)];
Xs = double(full(random(X,1,1000)));

%%
pc = POLYCHAOS(2,3);

figure(12)
clf
subplot(1,2,1)
multipdfsampleplot(Xs,'npts',30)
ax0 = axis;
cax0=caxis;
pause(2)

L=inf;
scanL=[];
scanXpc = cell(1,0);
%%
for i=1:10
    [Xpc0,L0] = PCidentification(pc,Xs,'nbestim',1e3,'rs',15,'fmin',10,'parametrized');
    if L0<L
        Xpc = Xpc0;
        L=L0;
        scanL=[scanL;L];
        scanXpc=[scanXpc , {Xpc}];
        scanL
        figure(12)
        subplot(1,2,2)
        multipdfplot(Xpc,'npts',30,'nbs',1e5)
        axis(ax0)
        caxis(cax0)
        pause(1)
    end
end



%% 

X = RVLOGNORMAL(1,0.3,0,'stat');
Xs = random(X,1,300);

Xpc = eicdfpcproject(Xs',8);

figure(100)
clf
pdfsampleplot(Xs,'r')
hold on

pdfplot(Xpc,'b','nbs',1e5,'ksdensity');

%%
X = RVNORMAL(0,0.3);
Xs = random(X,1,300);
X = RVNORMAL(1,0.3);
Xs = [Xs,random(X,1,600)];

Xpc = eicdfpcproject(Xs',20,[],'nbgauss',40);

figure(100)
clf
pdfsampleplot(Xs,'r')
hold on

pdfplot(Xpc,'b','nbs',1e5,'ksdensity');



