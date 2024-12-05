function P = calc_Pnode(elem,xnode)
% function P = calc_Pnode(elem,xnode)

R = calc_R(elem);

switch getindim(elem)
    case 2
        P = [R,zeros(2,1);zeros(1,2),1];
    case 3
        P = [R,zeros(3,3);zeros(3,3),R] ;
end
