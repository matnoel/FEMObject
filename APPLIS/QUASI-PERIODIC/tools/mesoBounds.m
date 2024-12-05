function bounds = mesoBounds(t,cellNum,iBndry,dist)
% bounds = mesoBounds(t,cellNum,iBndry,dist)
% _iBndry is index of boundary nodes in cellModel, as obtained from
%  getnumber(getnode(create_boundary(getCellModel(model))))
% _dist is optional
% _bounds is an array of size cellNumber*4 where each row i is
%  [inf sup inf_on_boundary sup_on_boundary] for cell i.
  
if nargin < 4
    dist = [] ;
end

order = t.order ;
cellNb = prod(cellNum) ;

% Get mesoscopic coordinates of one cell per phase.
if isempty(dist) % assume each cell is a different phase
    dist = num2cell(1:cellNb) ;
else
    dist(cellfun(@isempty,dist)) = [] ;
end
mesoCoord = cellfun(@(x) x(1),dist(:)) ;
mesoCoord = formatIndex(order,cellNum,mesoCoord) ;

% Evaluate at coordinates and find eigenvalues at every node
if iscell(t)
    assert(norm(t{1,2})==0 && norm(t{2,1})==0,'Non-diagonal not implemented')
    values1 = evalAtIndices(t{1,1},mesoCoord,1:order-1) ;
    values1 = double(values1)' ;
    values2 = evalAtIndices(t{2,2},mesoCoord,1:order-1) ;
    values2 = double(values2)' ;
    values = cat(3,values1,values2) ;
    eigenvalues = max(values,[],3) ;
elseif isa(t,'TuckerLikeTensor')
    if isa(t.space,'TSpaceVectors')
        eigenvalues = evalAtIndices(t,mesoCoord,1:order-1) ;
        eigenvalues = double(eigenvalues)' ;
    elseif isa(t.space,'TSpaceOperators')
        eigenvalues = zeros(t.sz(1,end)/2,prod(t.sz(1,1:end-1)),2) ;
        for i = 1:size(mesoCoord,1)
            ti = evalOperatorAtIndices(t,mesoCoord(i,:),mesoCoord(i,:),1:order-1) ;
            for k = 1:ti.sz(end)/2
                tik = evalOperatorAtIndices(ti,[2*k-1;2*k],[2*k-1;2*k],ti.order) ;
                tik = double(squeeze(full(tik),1:2)) ;
                eigenvalues(k,i,:) = sort(eig(tik)) ;
            end
        end
    end
end

% Get extreme eigenvalues across all nodes, for every cell
minVal = [min(eigenvalues(:,:,1),[],1) ; min(eigenvalues(iBndry,:,1),[],1)] ;
maxVal = [max(eigenvalues(:,:,end),[],1) ; max(eigenvalues(iBndry,:,end),[],1)] ;
        
% Store extrema at the right places
bounds = zeros(cellNb,4) ;
for i = 1:numel(dist)
    bounds(dist{i},1) = minVal(1,i) ;
    bounds(dist{i},2) = maxVal(1,i) ;
    bounds(dist{i},3) = minVal(2,i) ;
    bounds(dist{i},4) = maxVal(2,i) ;
end
end