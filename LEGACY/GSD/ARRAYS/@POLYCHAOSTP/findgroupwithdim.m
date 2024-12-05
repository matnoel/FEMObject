function [g,r]=findgroupwithdim(PC,i)
%function [g,r]=findgroupwithdim(PC,i)

g=[];r=[];
for k=1:getnbgroups(PC)
   
    [o,I] = ismember(i,PC.groups{k});
    if o
        g=k;
        r=I;
        break
    end
end
