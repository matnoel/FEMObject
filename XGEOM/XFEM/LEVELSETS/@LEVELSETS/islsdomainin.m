function [r,I]= islsdomainin(ls)
% function [r,i]= islsdomainin(ls)
r=false;I=[];
for i=1:length(ls.LS)
switch getnature(ls.LS{i})
    case 'domain'    
r= true;
I=[I,i];
end
end