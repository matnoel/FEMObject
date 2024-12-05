function A=purge(A,dim)
%      function A=purge(A,dim)
% Elimine les rangs = 0 en testant les dimensions dim
% Si dim non specifie, toutes les dimensions sont testees.
if nargin==1
    dim=1:A.dim;
end

keeprank=find(prod(cellfun(@(f) norm(f(:)),A.F(:,dim)),2)~=0);
A.F=A.F(keeprank,:);
A.alpha=A.alpha(keeprank);
A.m=length(A.alpha);