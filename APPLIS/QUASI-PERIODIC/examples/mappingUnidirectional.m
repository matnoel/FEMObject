approxTol = 1e-2 ;
barWidth = .5 ;
cellNb = 5 ;

model = QPModel('order',2,...
    'cellNum',[cellNb 1],...
    'cellSize',[1 5],...
    'elementSize',[.04 .04],...
    'tolSVD',1e-6,...
    'verbose',true) ;
tr = Truncator('tolerance',getTolSVD(model)) ;
patterns = struct('name',{'uniform','bar'},...
    'value',{1 99},...
    'size',{[] [barWidth 1]} ,...
    'center',{[] [.5 .5]},...
    'offset',{[] []}) ;
patternsTable = [1 1 ; 0 1] ;
cellCoord = getCellCoord(model) ;
microPatterns = drawCellPattern(cellCoord,patterns) ;
dist = {(1:cellNb)' ; setdiff((1:cellNb)',2)} ;
Kper = distributeMicro(model,dist,microPatterns*patternsTable) ;

%% Mapping

% Microscopic (patterns)
mapGPatterns = struct('name',{'uniform','bar','bar'},...
    'value',{1 1 1},...
    'size',{[] [.25 1] [.25 1]} ,...
    'center',{[] [.125 .5] [.875 1]},...
    'offset',{[] [] []}) ;
% mapGMicro = drawCellPattern(cellCoord(model),mapGPatterns) ;
locAll = microPatterns(:,1) ;
locBar = microPatterns(:,2)>0 ;
locBeforeBar = ~locBar & cellCoord(:,1)<.5 ;
locAfterBar = ~locBar & cellCoord(:,1)>.5 ;
mapGMicro = zeros(getNbCellDoF(model)*2,3) ;
mapGMicro(:,1) = 1 ;
mapGMicro(1:2:end,2:3) = [locBeforeBar locAfterBar] ;
mapG = distributeMicro(model,dist([1 1 1]),mapGMicro) ;
% mapMicro = [locAll cellCoord(:,1).*[locBeforeBar locBar locAfterBar]] ;
mapMicro = ones(getNbCellDoF(model)*2,4) ;
mapMicro(2:2:end,:) = 0 ;
mapMicro1 = ones(getNbCellDoF(model),1) ;
mapMicro2 = ones(getNbCellDoF(model),1) ;
mapMicro3 = ones(getNbCellDoF(model),1) ;
mapMicro2(locBeforeBar,1) = 0 ;
mapMicro2(locBar,1) = cellCoord(locBar,1) - min(cellCoord(locBar)) ;
mapMicro2(locAfterBar,1) = max(mapMicro2(locBar,1)) ;
locBeforeBar = mapMicro2(:,1)==0 ; % redefine to have
locAfterBar = mapMicro2(:,1)==max(mapMicro2(:,1)) ; % junction
mapMicro1(locBeforeBar,1) = cellCoord(locBeforeBar,1) ;
mapMicro1(~locBeforeBar,1) = max(mapMicro1(locBeforeBar,1)) ;
mapMicro3(locAfterBar,1) = cellCoord(locAfterBar,1) - min(cellCoord(locAfterBar)) ;
mapMicro3(~locAfterBar,1) = 0 ;
mapMicro(1:2:end,2:4) = [mapMicro1 mapMicro2 mapMicro3] ;
map = getCoord(model) ;
map = map{2} ;
map.space.spaces{end} = kron(map.space.spaces{end},[0;1]) ;
map = distributeMicro(model,dist([1 1 1 1]),mapMicro) + updateAllProperties(map) ;

% Mesoscopic (dilatations)
dilatation = 1.5*rand(cellNb+1,1)+.5 ; % U([0.5,2])
mapMeso = @(d) [ (0 + barWidth*((1:cellNb)'-1) + ... % as function of 
    (1-barWidth)*( cumsum([0;d(1:end-2)]) + cumsum([0;d(2:end-1)]))/2) ...% dilatations
    d(1:end-1) ones(cellNb,1) d(2:end)] ;
% invMap = map ;
map.space.spaces{1}(:,1:4) = mapMeso(dilatation) ;
% invMap.space.spaces{1} = mapMeso((dilatation).^(-1)) ;
% invMapG = mapG ;
mapG.space.spaces{1}(:,2:3) = [dilatation(1:end-1) dilatation(2:end)]-1 ;
% invMapG.space.spaces{1}(:,2:end) = (1+mapG.space.spaces{1}(:,2:end)).^(-1)-1 ;
detMapG = mapG ;
detMapG.space.spaces{end} = detMapG.space.spaces{end}(1:2:end,:) ;
detMapG = updateAllProperties(detMapG) ;
opDetMapG = detMapG ;
opDetMapG.space.spaces{end} = kron(opDetMapG.space.spaces{end},[1;1]) ;
opDetMapG = updateAllProperties(opDetMapG) ;
opDetMapG = toOperator(opDetMapG) ;
mapG = toOperator(mapG) ;

%% Periodic

% Assembling
pb1 = QPProblem(model,Kper,'corrector1',QPBC(4),'penalty',1e6) ;
[sol1,out1]=solve(pb1);

%% Matrix K2

K2 = Kper ;
K2.space.spaces{end} = kron(K2.space.spaces{end},[1;1]) ;
K2.space = updateProperties(K2.space) ;
K2 = toOperator(K2) ;
K2 = opDetMapG*(mapG*K2*mapG') ;
K2 = tr.truncate(K2) ;
pb2 = QPProblem(model,K2,'corrector1',QPBC(4),'penalty',1e6) ;
[sol2,out2]=solve(pb2);

%% Post-processing

figure
plot(model,sol1)
title('u_{per}')

figure
plot(model,detMapG)
title('\nabla \psi \cdot e_1')
figure
map1 = map ;
map1.space.spaces{end} = map1.space.spaces{end}(1:2:end,:) ;
map1 = updateAllProperties(map1) ;
plot(model,map1)
title('\psi \cdot e_1')

figure
coord = getDomainCoord(model) ;
mapCoord = [doubleQP(map1) coord(:,2)] ;
plot(model,Kper,'coord',mapCoord)

figure
plot(model,sol2)
title('u\circ\psi')

figure
plot(model,sol2,'coord',mapCoord)
title('u')

%% PPP

figure
colour = {'k','b','r','g','m','c'} ;
i = coord(:,2)==0 ;
hold on
for r = 1:map1.space.dim(1)
    y=kron(map1.space.spaces{1}(:,r),map1.space.spaces{2}(:,r));
    plot(coord(i,1),y(i),colour{r})
end
y = doubleQP(map1) ;
plot(coord(i,1),y(i),colour{r+1})
hold off
legend({'all','before bar','bar','after bar','total'})

figure
colour = {'k','b','r','g','m','c'} ;
i = cellCoord(:,2)==0 ;
hold on
for r = 1:map1.space.dim(1)
    y=map1.space.spaces{2}(i,r);
    plot(cellCoord(i,1),y,colour{r})
end
y = sum(map1.space.spaces{2}(i,:),2) ;
plot(coord(i,1),y,colour{r+1})
hold off
legend({'all','before bar','bar','after bar','total'})