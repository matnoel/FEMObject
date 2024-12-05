function t = tensorize(model,f,tolerance)
% t = tensorize(model,f,tolerance)

if nargin < 3
    tolerance = getTolSVD(model) ;
end

%% Cell Recursion

if iscell(f)
   assert(isfloat(f{1}),...
       'QPModel.tensorize: second input is float or cell array of float.')
   t = cell(size(f)) ;
   for i = 1:numel(f)
       t{i} = tensorize(model,f{i}) ;
   end
   return
end

%% Method

% For speed and clarity
tOrder = getOrder(model) ;
cellNb = getCellNb(model) ;
isVector = size(f,2) == 1 ;


% Build mesoscale indicators
if isVector
    mesoDoFNb = cellNb ;
    coordMeso = (1:mesoDoFNb)' ;
else % is operator
    mesoDoFNb = cellNb^2 ;
    [coordMeso(:,1),coordMeso(:,2)] = ind2sub(cellNb*[1 1],(1:mesoDoFNb)') ;
end
% Store each cell index (cells indices if operator) as a cell in array
coordMeso = mat2cell(coordMeso,ones(1,size(coordMeso,1)),size(coordMeso,2)) ;
mesoInd = mesoIndicator(model,coordMeso,~isVector) ;

% Build core
if tOrder == 2
    core = DiagonalTensor(ones(mesoDoFNb,1),tOrder) ;% each indicator is elementary
else % order 3
    indRank = cellfun(@(x) size(x,2),mesoInd(1,:)) ;
    totalRank = sum(indRank) ; % dimension of both meso subspaces
    core = zeros(totalRank,totalRank,mesoDoFNb) ; 
    for i = 1:mesoDoFNb
        currentRank = 1+sum(indRank(1:i-1)) ;
        range = currentRank:(currentRank+indRank(i)-1) ;
        core(range,range,i) = eye(indRank(i)) ;
    end
    core = FullTensor(core,3,[totalRank,totalRank,mesoDoFNb]) ;
end

% Build space
space = cell(tOrder,1) ;
% Mesoscale indicators
for o = 1:tOrder-1
    % Store each order in separate cell in cell array,
    % toward TSpaceOperators
    space{o} = [mesoInd{o,:}] ;
    if ~isVector ; space{o} = space{o}' ; end % as column for operators
end
% Microscale values
f = sunder(model,f) ; % duplicate node values on cells faces
if isVector
    % Fields over reference cell
    space{end} = reshape(f,[],cellNb) ;
    space = TSpaceVectors(space) ;
else % is operator
    % Split matrix into cell array, per mesoscale coord
    operatorSizes = diag(size(f)/cellNb)*ones(2,cellNb) ;
    cellOperators = mat2cell(f,operatorSizes(1,:),operatorSizes(2,:)) ;
    % Get microscale fields matching every mesoscale coord
    space{end} = cellOperators(:) ;
    space = TSpaceOperators(space) ;
end

% Build tensor
t = TuckerLikeTensor(core,space) ;

%% Compression

if tolerance > 0
    % Free memory
    clear core space cellOperators f mesoInd coordMeso
    
    % Compress
    tr = Truncator('tolerance',tolerance,'maxrank',mesoDoFNb) ;
    t = tr.truncate(t) ;
    
    % Ensure class (avoid CanonicalTensor)
    t = TuckerLikeTensor(t.core,t.space) ;
end

end

%% Previous method (no transfer index provided)

% function [t,fe2ts,ts2fe] = tensorize(model,f,coordOrIndex)
% % [t,fe2ts,ts2fe] = tensorize(model,f,coordOrIndex)
% 
% %TODO: use coord2indicator to get all indicators in one go ?
% 
% % For speed and clarity
% tOrder = getOrder(model) ;
% cellNb = getCellNb(model) ;
% cellNum = getCellNum(model) ;
% fColNb = size(f,2) ;
% 
% space = cell(tOrder,1) ;
% indexProvided = size(coordOrIndex,2) == 1 ;
% fe2ts = getFE2TS(model) ;
% if indexProvided || ~isempty(fe2ts)
%     %% Faster method: "FEM to tensor" index available
%     
%     % Fields over reference cell
%     if indexProvided % use provided index
%         fe2ts = coordOrIndex ;
%     end
%     f = f(fe2ts,:) ;
%     dofNb = size(f,1)/cellNb ;
%     if fColNb == 1 % f is scalar-valued
%         space{2} = reshape(f,dofNb,cellNb) ;
%     else % f is either vector-valued or an operator
%         space{2} = mat2cell(f,dofNb*ones(cellNb,1),fColNb) ;
%     end
%     
%     % Build mesoscale indicators
%     if tOrder==2 % order2
%         for i = 1:cellNb
%             mesoInd = coord2indicator(model,i) ;
%             if fColNb == 1 % f is scalar-valued
%                 space{1} = [space{1} mesoInd(:)] ;
%             elseif fColNb == getDim(model) % f is vector-valued
%                 space{1} = [space{1} ; {mesoInd(:)}] ;
%             else % f is an operator
%                 space{1} = [space{1} ; {diag(mesoInd(:))}] ;
%             end
%         end
%     else % order 3
%         for i = 1:cellNum(1)
%             for j = 1:cellNum(2)
%                 mesoInd = coord2indicator(model,[i j]) ;
%                 if fColNb == 1 % f is scalar-valued
%                     space{1} = [space{1} mesoInd(:,1)] ;
%                     space{2} = [space{2} mesoInd(:,2)] ;
%                 elseif fColNb == getDim(model) % f is vector-valued
%                     space{1} = [space{1} ; {mesoInd(:,1)}] ;
%                     space{2} = [space{2} ; {mesoInd(:,2)}] ;
%                 else % f is an operator
%                     space{1} = [space{1} ; {diag(mesoInd(:,1))}] ;
%                     space{2} = [space{2} ; {diag(mesoInd(:,2))}] ;
%                 end
%             end
%         end
%     end
% else
%     %% Slower method: no "FEM to tensor" available
%     domainCoord = coordOrIndex ;
%     mesoLower = mesoCoord(model,domainCoord,true) ;
%     mesoHigher = mesoCoord(model,domainCoord,false) ;
%     fe2ts = [] ;
%     
%     if tOrder==2 % order 2
%         for i = 1:cellNb
%             mesoInd = coord2indicator(model,i) ;
%             cellFE2TS = (mesoLower==i) | (mesoHigher==i) ;
%             fe2ts = [fe2ts ; cellFE2TS] ;
%             if fColNb == 1 % f is scalar-valued
%                 space{1} = [space{1} mesoInd(:)] ;
%                 space{2} = [space{2} f(cellFE2TS)] ;
%             elseif fColNb == getDim(model) % f is vector-valued
%                 space{1} = [space{1} ; {mesoInd(:)}] ;
%                 space{2} = [space{2} ; {f(cellFE2TS,:)}] ;
%             else % f is an operator
%                 space{1} = [space{1} ; {diag(mesoInd(:))}] ;
%                 space{2} = [space{2} ; {f(cellFE2TS,:)}] ;
%             end
%         end
%         
%     else % order 3
%         for i = 1:cellNum(1)
%             for j = 1:cellNum(2)
%                 mesoInd = coord2indicator(model,[i j]) ;
%                 cellFE2TS = all(mesoLower==[i j]) | all(mesoHigher==[i j]) ;
%                 fe2ts = [fe2ts ; cellFE2TS] ;
%                 if fColNb == 1 % f is scalar-valued
%                     space{1} = [space{1} mesoInd(:,1)] ;
%                     space{2} = [space{2} mesoInd(:,2)] ;
%                     space{3} = [space{3} f(cellFE2TS)] ;
%                 elseif fColNb == getDim(model) % f is vector-valued
%                     space{1} = [space{1} ; {mesoInd(:,1)}] ;
%                     space{2} = [space{2} ; {mesoInd(:,2)}] ;
%                     space{3} = [space{3} ; {f(cellFE2TS,:)}] ;
%                 else % f is an operator
%                     space{1} = [space{1} ; {diag(mesoInd(:,1))}] ;
%                     space{2} = [space{2} ; {diag(mesoInd(:,2))}] ;
%                     space{3} = [space{3} ; {f(cellFE2TS,:)}] ;
%                 end
%             end
%         end
%     end
% end
% 
% if nargout > 2 % because unique is costly
%     [~,ts2fe] = unique(fe2ts) ;
% end
% 
% if fColNb == 1
%     space = TSpaceVectors(space) ;
% else
%     space = TSpaceOperators(space) ;
% end
% core = DiagonalTensor(ones(cellNb,1),tOrder) ;
% t = TuckerLikeTensor(core,space) ;
% 
% % Compression
% tr = Truncator('tolerance',getTolSVD(model),'maxrank',cellNb) ;
% t = tr.truncate(t) ;
% 
% end