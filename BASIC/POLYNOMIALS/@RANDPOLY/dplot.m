function dplot(h,p,varargin)
% function dplot(h,p,varargin)
% plot of the derivative of a RANDPOLY h of degree p
% options of plot function 
%
% function plot(h,p,'xlim',xlim,'npts',n)
% x-axis limits for the plot
% n : number of points

xli = getcharin('xlim',varargin,getdomainborne(RANDVAR(h)));
n = getcharin('npts',varargin,200);

x = linspace(xli(1),xli(2),n);
hx = dpolyval(h,p,x);

varargin = delcharin('xlim',varargin);
varargin = delcharin('npts',varargin);

plot(x,hx,varargin{:});
xlim(xli);





