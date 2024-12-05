function pc = getpcgroups(pctp,i)

if nargin==1
pc = pctp.PCgroups;
else
pc = pctp.PCgroups(i);    
end


