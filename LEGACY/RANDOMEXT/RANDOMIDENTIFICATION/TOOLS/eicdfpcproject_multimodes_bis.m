function [Xpc,Ym] = eicdfpcproject_multimodes_bis(Xs,p,h,z,varargin)
% function Xpc = eicdfpcproject(Xs,p,h,z)
% p : ordre du chaos
% h: RANDPOLY
% z : vecteur permettant de separer les echantillons en differents blocks
% 
% See also eicdfpcproject, separatesamples

if nargin<=2 || isempty(h)
    h = POLYHERMITE();
end
m = length(z)+1;
[P,x,Xblock] = separatesamples(Xs,z,varargin{:});

PC = POLYCHAOS(RANDPOLYS(POLYFE(x),h),[0,p],'typebase',2);

ng = getcharin('nbgauss',varargin);
if length(ng)<2
    ng = [1,ng];
end
Xpc = decompfun(PC,ng,[],@(xi) myfun(xi,RANDVAR(h),Xblock,x));

ind = getindices(PC);
Ym = cell(1,m);
pc = POLYCHAOS(h,p);

for i=1:m
    rep = find(ind(:,1)==i-1);
    Ym{i} = double(getpccompo(Xpc,rep));
    Ym{i} = PCMATRIX(Ym{i},[1,1],pc);
end

function y = myfun(xi,rv,Xblock,x,P)

xi2 = cdf(rv,xi(:,2));
xi1 = xi(:,1);
y = zeros(size(xi,1),1);

for i=1:size(xi,1)
    k =  max(find(x<xi1(i))); 
    if isempty(k)
        k=1;
    end
    if k==1
        rep = find(xi1(i)<=x(2));
    else
        rep = find(xi1(i)>x(k) | xi1(i)<=x(k+1));
    end
    y(i) = eicdf(Xblock{k},xi2(i));  
end

return

