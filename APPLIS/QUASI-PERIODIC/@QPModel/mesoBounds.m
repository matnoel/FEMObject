function bounds = mesoBounds(model,t,dist)
% bounds = mesoBounds(model,t,dist)
% _dist is optional
% _bounds is an array of size cellNumber*4 where each row i is
%  [inf sup inf_on_boundary sup_on_boundary] for cell i.

if nargin < 3
    dist = [] ;
end

iBndry = getnumber(getnode(create_boundary(getCellModel(model)))) ;
bounds = mesoBounds(t,getCellNum(model),iBndry,dist) ;

end