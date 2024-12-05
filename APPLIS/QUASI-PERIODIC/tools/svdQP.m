function y = svdQP(m,cellNum,tol)
% y = svdQP(m,cellNum,tol)
% Convert vector of V(I) (order 2 format) into vector of V(I)\otimes V(J)
% (order 3 format). Alternatively, converts operator over V(I) into 
% operator over V(I)\otimes V(J). SVD compression is used.
% Cell recursion is implemented.
%
% m : column of length prod(cellNum) or square matrix of size prod(cellNum).
% cellNum : [n1 n2] with ni number of cells along dimension i.
% tol : [optional] tolerance for truncated SVD. Default 0 (no truncation).
% y : vectors of (or operators over) V(I) and V(Y), in cells as in 
% TSpaceVectors (TSpaceOperators) format. Cell array of such if m is a cell
% array.

if nargin < 3
    tol = 0;
end

%% Cell recursion

if iscell(m)
    y = cell(size(m)) ;
    for n = 1:numel(m)
        y{n} = svdQP(m,cellNum,tol) ;
    end
end

%% Method

if issparse(m) ; m = full(m) ; end % safety for reshape

isOperator = size(m,2)>1 ;
if isscalar(m) 
    warning('Degenerate case: size 1x1, assumed not to be operator') 
end

if isOperator
    % Re-organize M_{\lambda\lambda'} into M_{ii'jj'} as 2D matrix
    m = reshape(m,[cellNum cellNum]) ; % here is M_{iji'j'}
    m = permute(m,[1 3 2 4]) ; % here is M_{ii'jj'}
    m = reshape(m,cellNum.^2) ; % here is the 2D matrix
else
    m = reshape(m,cellNum) ;
end

[u,s,v] = svd(m,'econ') ;

% Compute SVD truncation errors and select rank according to tolerance
if verLessThan('matlab','8.2') % compatibility (<R2013b)
    truncErr = sqrt(flipdim(cumsum(flipdim(s.^2,1)),1)/sum(s.^2));
else
    truncErr = sqrt(flip(cumsum(flip(s.^2)))/sum(s.^2));
end
truncErr = [truncErr(2:end) ; 0] ; % shift "ranks errors" from 0:n-1 to 1:n
rank = find(truncErr<tol,1) ; % take lowest index satisfying tolerance

% Store in TSpaceVectors format
% Every vector is a column of an array. Each order has its array stored
% in a separate cell (TSpaceVectors cell array format)
y = {zeros(size(u,1),rank) ; zeros(size(v,1),rank)} ;
for i = 1:rank
    y{1}(:,i) = u(:,i)*s(i,i).^0.5 ;
    y{2}(:,i) = v(:,i)*s(i,i).^0.5 ;
    % Output expected positive so, in case both u and v are negatives:
    if all(y{1}(:,i)<=0) && all(y{2}(:,i)<=0)
        y{1}(:,i) = -y{1}(:,i) ;
        y{2}(:,i) = -y{2}(:,i) ;
    end
end

% Output expected positive so, in case both u and v are negatives:
% if all(cellfun(@(x) any(any(x<0)),y)) % unless u and v are all positive
%     bothNegatives = find(all(y{1}<=0,1) & all(y{2}<=0)) ;
%     y{1}(:,bothNegatives) = abs(y{1}(:,bothNegatives)) ;
%     y{2}(:,bothNegatives) = abs(y{2}(:,bothNegatives)) ;
% end


% Store in TSpaceOperators format
% For each order, reshape every column into matrix and store in
% separate cell (TSpaceOperators cell array format)
if isOperator
    yy = cell(2,1) ;
    for c = 1:size(y{1},2)
        yy{1} = [yy{1} ; {reshape(y{1}(:,c),cellNum(1)*[1 1])}] ;
        yy{2} = [yy{2} ; {reshape(y{2}(:,c),cellNum(2)*[1 1])}] ;
    end
    y = yy ;
end
end

