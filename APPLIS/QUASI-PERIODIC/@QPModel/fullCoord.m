function x = fullCoord(model,i,y)
% x = fullCoord(model,i,y)

if size(i,2) > 2 && size(i,1) == 1 % is a row of indices
    i = i' ;
end

if size(i,2) == 1 % is a column of indices
    cellNum = getCellNum(model) ;
    i = [1+mod(i-1,cellNum(1)) ceil(i/cellNum(1))] ;
end

x = y + (i-1)*diag(getCellSize(model)) ;

end