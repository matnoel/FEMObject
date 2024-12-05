function level = getlevelset(level,k)
% function level = getlevelset(level,k)

if nargin==2 &&  k~=level.number
   error('wrong levelset number') 
end