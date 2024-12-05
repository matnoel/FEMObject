function [r,I]= islscrackin(ls)
% function [r,i]= islscrackin(ls)
r=false;I=[];
for i=1:length(ls.LS)
switch getnature(ls.LS{i})
    case 'crack'    
r= true;
I=[I,i];
end
end