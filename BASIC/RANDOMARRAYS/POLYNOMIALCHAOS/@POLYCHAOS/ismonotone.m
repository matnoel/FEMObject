function r=ismonotone(PC)
% function r=ismonotone(PC)
% 
% Check whether or not the multi-index set of POLYCHAOS PC is monotone (or downward closed or lower)
%
% See also POLYCHAOS/isdownwardclosed, POLYCHAOS/islower

ind=PC.indices;
M=PC.M;
P=length(PC);
h=PC.RANDPOLYS;

r = true;
for i=1:P
    p = ind(i,1:M);
    PC_rect = POLYCHAOS(h,p,'typebase',2);
    ind_rect = getindices(PC_rect);
    [ok,rep]=ismember(ind_rect(:,1:M),ind(:,1:M),'rows');
    if ~all(ok)
        r = false;
        return
    end
end

end