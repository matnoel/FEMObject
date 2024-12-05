function c = LSCRACK(ls1,ls2,ls3)
% function c = LSCRACK(ls1,ls2,ls3)

if nargin==0
    c = struct();
    ls = cell(0,1);
    c.number=[];
    c = class(c,'LSCRACK',LEVELSETS(ls));
else
    c = struct();
    ls=cell(1,nargin);
    ls{1} = ls1;
    if nargin>=2
        ls{2} = ls2;
    end
    if nargin>= 3
        ls{3} = ls3;
    end
    c.number = [];
    c = class(c,'LSCRACK',LEVELSETS(ls));
end
