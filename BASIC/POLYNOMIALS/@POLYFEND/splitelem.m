function H = splitelem(H,E,j)
% function H = splitelem(H,E)
% decoupage des elements E
%
% function H = splitelem(H)
% decoupage de tous les elements
% 
% function H = splitelem(H,E,j)
% decoupage des elements E selon las dimension  j
%
% function H = splitelem(H,[],j)
% decoupage de tous les elements selon la dimension  j

if nargin==1 || isempty(E)
    E = 1:getnbelem(H);
end
M = getM(H);
Z = H.e;
for k=1:length(E)
    e = E(k);
    Ze = Z{e};

    if nargin<=2
        Znew = cell(1,2^M);
        waynew = calc_multiindices(M,2);
        for l=1:2^M
         Znew{l}.way = [Ze.way;waynew(l,:)];  
         Znew{l}.order = Ze.order+1;   
        end
    else
        Znew = cell(1,2);  
        waynew = ones(2,M);
        waynew(2,j) = 2;      
        if numel(Ze.order)==1
            Ze.order = repmat((0:Ze.order)',1,size(Ze.way,2));
        end
        ordernew = Ze.order(end,:);
        ordernew(j) = ordernew(j)+1;
        for l=1:2
        Znew{l}.way = [Ze.way;waynew(l,:)];  
        Znew{l}.order = [Ze.order;ordernew];
        Znew{l}.order(j) = Znew{l}.order(j)+1;  
        end
        
    end
    Z = [Z(1:e-1),Znew,Z(e+1:end)];
end

%H=Z;
H = POLYFEND(M,Z);


