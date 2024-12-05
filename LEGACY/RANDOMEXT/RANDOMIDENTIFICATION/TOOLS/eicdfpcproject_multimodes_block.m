function [Xpc,Ym] = eicdfpcproject_multimodes_block(Xblock,P,p,h,varargin)
% function Xpc = eicdfpcproject(Xblock,p,h,varargin)
% Xblock : cell contenant les blocks de realisations
% P : vecteur donnant les probabilites de chaque evenement
% p : ordre du chaos
% h: RANDPOLY
% eicdfpcproject(p,h,Xblock,P,'nbgauss',n)
% n : nb de points de gauss pour la quadrature

N = 0;
m = length(Xblock);
%P=zeros(1,m);
for i=1:m
    N = N + numel(Xblock{i});
%    P(i) = numel(Xblock{i});
end
%P = P/N;

if nargin<=2 || isempty(h)
    h = POLYHERMITE();
end

x = zeros(1,length(P)+1);
for i=1:length(P)
    x(i+1) = x(i)+P(i);
end 

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






