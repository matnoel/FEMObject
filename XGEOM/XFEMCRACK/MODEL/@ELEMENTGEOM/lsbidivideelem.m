function [elemin,nodeplus,xnodein] = lsbidivideelem(elem,ls1,ls2,node)


[elemin,elemout,nodeplus,xnodein,xnodeout] = newlsdivideelem(elem,ls1,node);

elemin = [elemin,elemout];
xnodein = [xnodein,xnodeout];

keyboard
%for k=1:length(elemin)
%calc_N(elemin{k},'nbddlpernode',1)