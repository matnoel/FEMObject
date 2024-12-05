function s=numel(ls,varargin)
% function s=numel(ls)
% s est le nombre de levelset
if nargin==1
if isa(ls.value,'MULTIMATRIX')
s=length(ls.value);
elseif isa(ls.value,'double')
s=size(ls.value,2);
end

else
s=1; % pour quand on fait subsref avec {}
end