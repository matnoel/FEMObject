function [B,detJ,DN] = calc_gradient(elem,xnode,xgauss)
% function [B,detJ,DN] = calc_gradient(elem,xnode,xgauss)

nbelem=getnbelem(elem);

[DN,detJ] = calc_DN(elem,xnode,xgauss);

if getnbddlpernode(elem)==1
    nbddl = getnbnode(elem);
    B=DN;
elseif getnbddlpernode(elem)>=getindim(elem)
    nbddl = getnbnode(elem)*getindim(elem);
    switch getindim(elem)
        case 1
            B=DN;
        case 2
            B=zerosND([4,nbddl,sizeND(DN)]);
            B(1,1:2:end)=DN(1,:);
            B(2,2:2:end)=DN(1,:);
            B(3,1:2:end)=DN(2,:);
            B(4,2:2:end)=DN(2,:);

        case 3
            B=zerosND([9,nbddl,sizeND(DN)]);
            B(1,1:3:end)=DN(1,:);
            B(2,2:3:end)=DN(1,:);
            B(3,3:3:end)=DN(1,:);
            B(4,1:3:end)=DN(2,:);
            B(5,2:3:end)=DN(2,:);
            B(6,3:3:end)=DN(2,:);
            B(7,1:3:end)=DN(3,:);
            B(8,2:3:end)=DN(3,:);
            B(9,3:3:end)=DN(3,:);

    end
else
    error('calc_gradient non defini pour ce type de ddl aux noeuds')
end
