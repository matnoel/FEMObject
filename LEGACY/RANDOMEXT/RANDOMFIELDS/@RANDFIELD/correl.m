function co=correl(RF,x,y)
if nargin==1
 co = RF.correl;
elseif nargin==3
co=eval(RF.correl,x,y);
end
