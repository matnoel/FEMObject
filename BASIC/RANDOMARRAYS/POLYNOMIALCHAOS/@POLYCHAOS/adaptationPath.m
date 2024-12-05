function apath = adaptationPath(PC)
%function apath = adaptationPath(PC)
% creates a boolean matrix apath with column k 
% giving the basis functions with multiindices 
% with sum less than k-1 

 apath = false(PC.P+1,max(PC.n)); 
 for k=1:max(PC.n)
 apath(PC.indices(:,end)<=k-1,k)=true;
 end
