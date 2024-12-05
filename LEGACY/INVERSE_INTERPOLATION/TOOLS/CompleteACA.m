function [PivotValue,PivotElement,Rows,Cols,ifail,INFNORM] = CompleteACA(A,tol)
% Adaptive Cross Approximation with complete pivoting. This command is 
% completely analogous to Gaussian elimination with complete pivoting. We
% adaptively find the rank of the approximant in this command. 

% Copyright 2013 by The University of Oxford and The Chebfun Developers. 
% See http://www.maths.ox.ac.uk/chebfun/ for Chebfun information.

if nargin<2
    tol=0;
end

% Set up output variables.
[nx,ny]=size(A);
width = min(nx,ny);         % Use to tell us how many pivots we can take.
PivotValue = zeros(1);      % Store an unknown number of Pivot values.
PivotElement = zeros(1,2);  % Store (j,k) entries of pivot location.
ifail = 1;                  % Assume we fail.
factor = 4*(tol>0);         % ratio between size of matrix and no. pivots. If tol = 0, then do full no. of steps.

% Main algorithm
zrows = 0;                  % count number of zero cols/rows.
[ infnorm , ind ]=max( abs ( A(:) ) );
[ row , col ]=myind2sub( size(A) , ind);
scl = infnorm;
INFNORM=infnorm;


while ( infnorm > tol*scl )   && ( zrows < width ) %% && ( zrows < width / factor) )
    tmprow = A( row , : ) ;
    Rows(zrows+1,:)=tmprow;
    tmpcol = A( : , col ) ;    % Extract the columns out
    Cols(:,zrows+1) = tmpcol;
    PivVal = A(row,col);
    A = A - tmpcol*(tmprow./PivVal);    % Rank one update.
    
    % Keep track of progress.
    zrows = zrows + 1;                  % One more row is zero.
    PivotValue(zrows) = PivVal;         % Store pivot value.
    PivotElement(zrows,:)=[row col];    % Store pivot location.
    
    %Next pivot.
    [ infnorm , ind ]=max( abs ( A(:) ) );  
    [ row , col ]=myind2sub( size(A) , ind );
    INFNORM = [INFNORM infnorm];
end

if infnorm <= tol*scl, ifail = 0; end  % We didn't fail.

if ifail == 0, return; end


end




function [row col] = myind2sub(siz,ndx)
% My version of ind2sub. In2sub is slow because it has a varargout. Since
% this is at the very inner part of the constructor and slowing things down
% we will make our own. 
% This version is about 1000 times faster than MATLAB ind2sub. 

vi = rem(ndx-1,siz(1)) + 1 ; 
col = (ndx - vi)/siz(1) + 1;
row = (vi-1) + 1;

end