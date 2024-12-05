function [Xpc,Ym,P] = eicdfpcproject_multimodes(Xs,p,h,z,varargin)
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

ng = getcharin('nbgauss',varargin);
if length(ng)>2
    error('nbgauss : rentrer un entier')
end
PC = POLYCHAOS(RANDPOLYS(POLYFE(x),h),[0,p],'typebase',2);
ind = getindices(PC); 
Xpc = zeros(1,length(PC));
Ym=cell(1,m);

for i=1:m
    Ym{i} = eicdfpcproject(Xblock{i},p,h,'nbgauss',ng);
    rep = find(ind(:,1)==i-1);
    Xpc(rep)=sqrt(P(i))*double(Ym{i});
end

Xpc = PCMATRIX(Xpc,[1,1],PC);






