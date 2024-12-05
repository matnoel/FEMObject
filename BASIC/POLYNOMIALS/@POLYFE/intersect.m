function H = intersect(H1,H2)
% function H = intersect(H1,H2)

if ~ismember(H1,H2)
    error('les dimensions stochastiques doivent correspondre')
end
if isa(H1,'POLYFE') && isa(H2,'POLYFE')
    param1 = getparam(H1);
    param2 = getparam(H2);
    
    I = sortrows([param1.I;param2.I],1);
    
    % I = sort([param1.I,param2.I]);
    I(find(I(2:end,1)-I(1:end-1,1)<eps)+1,:)=[];
    I = sortrows(I,2);
    I(find(I(2:end,2)-I(1:end-1,2)<eps)+1,:)=[];
    
    H = POLYFE(I,min(param1.p,param2.p),'elem');
    
elseif isa(H1,'POLYFE')
    H=H1;
else
    H=H2;
end