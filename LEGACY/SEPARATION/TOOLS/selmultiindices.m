function [I,Ind] = selmultiindices(siz,type,varargin)
%function [I,Ind] = selmultiindices(s,type,varargin)
% selection of multiindices in [1...s(1)]x[1...s(2)]x[1...s(end)]
% 
% I: boolean indexing
% Ind: multiindices which are kept
%
% function I = selmultiindices(dim,s,'diagonaldistance',l)
% distance to the diagonal elements
% example in 2D: an index (i,j) is selected if abs(i-j)<=l
%                 l=0 only select the diagonal terms
Ind = cell(1,length(siz));
[Ind{:}]= ind2sub(siz,1:prod(siz));
Ind = vertcat(Ind{:})';
I = Ind-repmat(min(Ind,[],2),1,size(Ind,2));
switch type
case 'diagonaldistance'
    level=varargin{1};
    I = sum(I,2);
%I = find(I<=prod(siz));
    I = I<=level;
otherwise
    error('not programmed')
end

fprintf('percentage of indices kept = %.0f \n',nnz(I)/numel(I)*100)
Ind = Ind(I,:);
Ind = mat2cell(Ind,size(Ind,1),ones(1,size(Ind,2)));
I = sub2ind(siz,Ind{:});
