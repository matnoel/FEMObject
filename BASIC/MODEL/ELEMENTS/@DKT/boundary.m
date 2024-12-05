function [segbord,nodebord,xnodeseg] = boundary(elem,node)
% function [segbord,nodebord,xnodeseg] = boundary(elem,node)

[segbord,nodebord,xnodeseg] = boundary(elem.TRI3,node);

segbord = SEG2DK(segbord);

segbord = update_lsdata(segbord,elem);
