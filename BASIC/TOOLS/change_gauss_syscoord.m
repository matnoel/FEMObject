function gaussParent = change_gauss_syscoord(gaussChild,elemChild,elemParent)
% function gaussParent = change_gauss_syscoord(gaussChild,elemChild,elemParent)
%   Change coordinate system of Gauss points structure gaussChild
% (associated with elemChild) into the local coordinate system of
% elemParent. 
%   The third dimension of gaussParent.coord and gaussParent.w
% is the number of elements in elemParent (completed with zeros).
%   This was tested and implemented for BILINFORMBOUNDARY/eval_elem, with 
% QUA4 as elemParent and SEG2 as elemChild, the latter created from the
% former with create_boundary(modelParent,'withparent') (and possibly
% intersection with edges).
% 
% Note of Gauss points storage organisation
% gauss.coord's size seems to be (this is guess work): 
% #dof/node x #spatialDimensions x #elements x #GaussPoint
% gauss.w is consistently organised as
% #dof/node x 1 x #elements x #GaussPoint
% Size along the third dimension is often 1 since the Gauss points are the 
% same in every element of the group (in local coordinate system). This is
% not conserved by this function.
% Problems may occur in case of only one gauss point per element, because
% MATLAB collapses trailing dimensions of size 1; e.g. zeros(2,2,1) returns
% a 2-by-2 matrix. This issue gave rise to the conditional structures
% below.

connecChild = getconnec(elemChild) ;
connecParent = getconnec(elemParent) ;
xGauss = gaussChild.coord ;
xRefChild = nodelocalcoord(elemChild) ;
xRefParent = nodelocalcoord(elemParent) ;

% Permute to copy gauss.coord organisation (as in permutegaussND)
xRefChild = MYDOUBLEND(permute(xRefChild,[4,2,1,3])) ;
xRefParent = MYDOUBLEND(permute(xRefParent,[4,2,1,3])) ;

hChild = xRefChild(:,:,2,1)-xRefChild(:,:,1,1) ;
xi = (xGauss-xRefChild(:,:,1,1))./hChild ;
[c2pElemList,c2pRefNodesList] = find_edges(connecChild,connecParent) ;
c2pElemList = [getnumelem(elemChild) c2pElemList] ;

% Store duplicates (nth children with n>1) separately
c2pElem = [] ;
c2pRefNodes = [] ;
dup = (1:size(c2pElemList,1))' ; % potential duplicates
while ~isempty(dup)
    [~,o2u] = unique(c2pElemList(dup,2),'stable') ; % get original to unique index
    c2pElem = [c2pElem ; {c2pElemList(dup(o2u),:)}] ;
    c2pRefNodes = [c2pRefNodes ; {c2pRefNodesList(dup(o2u),:)}] ;
    dup(o2u) = [] ;
end
% Read: c2pElem{n}(k,2) has c2pElem{n}(k,1) as its nth child (for any k).

gaussParent = gaussChild ;
szCoord = size(gaussParent.coord) ;
if size(szCoord,2)<4 % if, e.g., only one gauss point per element
    szCoord = [szCoord ones(1,4-size(szCoord,2))];
end
szCoord(2:4) = [2 size(connecParent,1) szCoord(4)*numel(c2pElem)] ;
% szCoord = [szCoord(1) 2 size(connecParent,1) szCoord(end)*numel(c2pElem)] ; % better ?
allCoord = MYDOUBLEND(zeros(szCoord)) ;
szW = size(gaussParent.w) ;
if size(szW,2)<4 % if, e.g., only one gauss point per element
    szW = [szW ones(1,4-size(szW,2))];
end
nGauss = szW(4) ;
szW(3:4) = [size(connecParent,1) szW(4)*numel(c2pElem)] ;
% szW = [szW(1:2) size(connecParent,1) szW(end)*numel(c2pElem)] ; % better?
allW = MYDOUBLEND(zeros(szW)) ;

for k = 1:numel(c2pElem)
    hParent = xRefParent(:,:,c2pRefNodes{k}(:,2),1) ...
        - xRefParent(:,:,c2pRefNodes{k}(:,1),1);
    % Adapt to number of Gauss points per element
    currentPoints = (nGauss*k-(nGauss-1)):(nGauss*k) ;
    allCoord(:,:,c2pElem{k}(:,2),currentPoints) = ...
        xRefParent(:,:,c2pRefNodes{k}(:,1),1) + xi*hParent ;
%         + xi(:,:,c2pElem{k}(:,1))*hParent ; % check subscripts
    allW(:,:,c2pElem{k}(:,2),currentPoints) = ...
        gaussChild.w*(norm(hParent)/double(norm(hChild))) ;
end

gaussParent.coord = allCoord ;
gaussParent.w = allW ;
gaussParent.nbgauss = size(gaussParent.coord,4) ;

end
%TODO: Handle several edges to assemble over whole boundary model. 
% Best would be to handle any edge orientation.