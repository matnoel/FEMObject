function [r,I]= islsmaterialin(ls)
% function [r,i]= islsmaterialin(ls)
r=false;I=[];
for i=1:length(ls.LS)
switch getnature(ls.LS{i})
    case 'material'    
r= true;
I=[I,i];
end
end