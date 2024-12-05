function [i,y] = mesoMicroCoord(model,x)
% [i,y] = mesoMicroCoord(model,x)
% Boundary behaviour: lower cell indices are chosen.

i = mesoCoord(model,x,false) ;
y = microCoord(model,x) ;

end

function i = mesoCoord(model,x,upperValue)
% i = mesoCoord(model,x,upperValue)

if nargin < 3
    upperValue = false ; % Choose upper value on inter-cell boundary
    % this is the behaviour consistent with microCoord, such that
    % fullCoord and mesoMicroCoord are reciprocal
end

cellNum = getCellNum(model) ;
x = x*diag(getCellSize(model).^-1) ; % process as relative "x/cellSize"

if upperValue % Pre-process to get higher value on inter-cell boundary
    rep = mod(x,1) == 0 ; % get boundary values subscripts
    x(rep) = x(rep)+1 ; % set to higher value
end

i = ceil(x) ; % return lower values on boundaries if no pre-processing
i = max(1,i) ; % set zero values to 1
i = [min(cellNum(1),i(:,1)) min(cellNum(2),i(:,2))] ; % set cellNb+1 values to cellNb

% Format for order 2
if getOrder(model) == 2
    i = i(:,1) + cellNum(1)*(i(:,2)-1) ;
end

end

function y = microCoord(model,x)
% y = microCoord(model,x)

cSize = getCellSize(model) ;
x = x*diag(cSize.^-1) ; % process as relative "x/cellSize"
y = x-fix(x) ; % get decimal part
iCorrect = y==0 & x~=0 ; % only original zeros are accepted
y(iCorrect) = 1 ; % others are set to one (boundary choice consistency)
y = y*diag(cSize) ; % revert first operation

end