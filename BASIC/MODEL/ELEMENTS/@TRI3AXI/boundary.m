function [segbord,nodebord,xnodeseg] = boundary(elem,node)
% function [segbord,nodebord,xnodeseg] = boundary(elem,node)

[segbord,nodebord,xnodeseg] = boundary(elem.TRI3,node);
segbord = convert(segbord,node,'SEG2AXI');
