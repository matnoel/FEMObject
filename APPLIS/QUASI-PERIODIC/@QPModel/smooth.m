function continuous = smooth(model,discontinuous)
% continuous = smooth(model,discontinuous)

% Safety
if length(discontinuous) == getNbDomainDoF(model)
    continuous = discontinuous ;
    return
elseif length(discontinuous) ~= getNbTotalDoF(model)
    error('Input has incorrect length.')
end

isOperator = size(discontinuous,1) == size(discontinuous,2) ;

d2c = getDiscon2Con(model) ;

if isOperator
    d2c(d2c>0) = 1 ; % use sum instead of average
    continuous = d2c*discontinuous*d2c' ;
else
    continuous = d2c*discontinuous ;
end

% % Legacy code
%
% index = (1:length(discontinuous))' ;
% duplicates = setdiff(index,d2c) ;
% coord = getDomainCoord(model) ;
% tol = 1e-9 ;
% 
% processed = zeros(numel(duplicates),1) ;
% for i = duplicates' %
%     if ~any(processed(duplicates==i))
%         coordDiff = [coord(:,1)-coord(i,1), coord(:,2)-coord(i,2)] ;
%         identicals = find(max(abs(coordDiff),[],2)<tol) ;
%         if isOperator % sum across rows then columns
%             smoothValues = sum(discontinuous(identicals,:),1) ;
%             discontinuous(identicals,:) = repmat(smoothValues,numel(identicals),1) ;
%             smoothValues = sum(discontinuous(:,identicals),2) ;
%             discontinuous(:,identicals) = repmat(smoothValues,1,numel(identicals)) ;
%         else % average across rows
%             smoothValues = mean(discontinuous(identicals,:)) ;
%             discontinuous(identicals,:) = repmat(smoothValues,numel(identicals),1) ;
%         end
%         % For speed
%         processed(duplicates==i) = 1 ; % keep track of processed duplicates
%     end
% end
% 
% continuous = discontinuous(d2c,:) ;
% if isOperator
%     continuous = continuous(:,d2c) ;
% end
end