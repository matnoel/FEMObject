function bx = condtransfer(a,b,ax,y)
% function bx = condtransfer(a,b,ax,y)
% a : RANDVAR
% b : CONDRANDVAR b depend de M variables aleatoires (b1,b2,...,bM)
% ax : n-by-1 double contenant les realisations de a
% y : n-by-M double contenant les realisations de (b1,b2,...,bM)
%
% bx  = F_(b|b1x,b2x,...bMx)^-1(F_(a)(ax))
%
% ------------------------------------
% function bx = condtransfer(a,b,ax)
% a : RANDVARS , M+1 variables aleatoires independantes a=(a1,a2,...aM+1)
% b : CONDRANDVAR b depend de M variables aleatoires (b1,b2,...,bM)
%
% ax : n-by-(M+1) double , ax(:,i) (note axi) sont les realisations de ai
%
% on a alors
% bx1 = F_(b1)^-1(F_(a1)(ax1))
% bx2 = F_(b2|bx1)^-1(F_(a2)(ax2))
% bx3 = F_(b3|bx1,bx2)^-1(F_(a3)(ax3))
% ...
% bx  = F_(b|bx1,bx2,...bxM)^-1(F_(aM+1)(axM+1))
%

if isa(a,'RANDVARS')
    if length(a)~=length(b.Y)+1
        error('les dimensions stochastiques doivent correspondre')
    end
    if isa(b.Y{1},'CONDRANDVAR')
        error('la premiere variable ne doit pas etre conditionnelle')
    end
    bx=zeros(size(ax));
    bx(:,1)=transfer(a{1},b.Y{1},ax(:,1));
    for i=2:length(b.Y)
        if isa(b.Y{i},'RANDVAR')
            bx(:,i)=transfer(a{i},b.Y{i},ax(:,i));
        elseif isa(b.Y{i},'CONDRANDVAR')
            bx(:,i)=unitransfer(a{i},b.Y{i},ax(:,i),bx(:,1:i-1));
        end
    end
    
    bx = unitransfer(a{end},b,ax(:,end),bx(:,1:end-1));
    
elseif isa(a,'RANDVAR')
    
    bx = unitransfer(a,b,ax,y);
    
end


function bx=unitransfer(a,b,ax,y)
% function bx = condtransfer(a,b,ax,y)
% a : RANDVAR
% b : CONDRANDVAR
%
% ax : n-by-1 double
% y : n-by-M double o� M est le nombre de de variables dont depend b
% on note y = (y1,y2,...,yM)
% bx = F_(b|y1,...,yM)^-1(F_(a)(ax))

for i=1:length(b.Y)
    b.Y{i} = y(:,i);
end
for i=1:length(b.funparam)
    if isa(b.funparam{i},'function_handle') | isa(b.funparam{i},'inline')
        b.funparam{i}=b.funparam{i}(b.Y{:});
    end
end
b = b.X(b.funparam{:});

bx=icdf(b,cdf(a,ax));

return