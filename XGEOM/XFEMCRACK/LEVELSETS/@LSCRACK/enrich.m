function c = enrich(c,typesupport,typetip,varargin)
% function ls = enrich(ls,typesupport,typetip,varargin)
% typesupport : 1 enrichissement saut classique
%               2 enrichissement localise (laisse leur signification intiale aux ddl FEM)
% typetip : enrichissement en pointe

if nargin==1 || isempty(typesupport)
   typesupport = 1; 
end
if nargin==2 || isempty(typetip)
   typetip = 0;
end

c.LEVELSETS{1}=enrich(c.LEVELSETS{1},typesupport);
for k=1:getnbtip(c)
c.LEVELSETS{1+k}=enrich(c.LEVELSETS{1+k},typetip);    
end

