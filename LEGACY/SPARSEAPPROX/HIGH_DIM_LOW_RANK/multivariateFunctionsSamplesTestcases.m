function [rs,fs,RV] = multivariateFunctionsSamplesTestcases(cas,N,varargin)
% function [rs,fs,RV] = multivariateFunctionsSamplesTestcases(cas ,N,varargin)
% evaluate N samples of a function f specified by k
% cas = 1: f(x) = x1 + ... + xd (varargin{1}=d), xi uniform variables on (-1,1)
% cas = 2: orange testcase
% cas = 3: oscillatory f(x) = cos(2piw+sum_i c_i x_i) (varargin = {d,w,c})
% cas =4 : gaussian
% cas = 5 : borehole 
switch cas
    case 1
d = varargin{1};        
RV = RANDVARS(RVUNIFORM(-1,1),d);
rs = random(RV,N,1);
rs = cell2mat(rs);
fs = sum(rs,2);
    case 2
[rs,fs,RV] = orange(N)  ;  
    case 3 % oscillatory
d = varargin{1};        
RV = RANDVARS(RVUNIFORM(-1,1),d);        
rs = random(RV,N,1);
rs = cell2mat(rs);
%fs = cos(2*pi*varargin{2}+sum(rs*diag(varargin{3}),2));  
fs = sin(sum(rs,2));  

    case 4 % gaussian
d = varargin{1};        
RV = RANDVARS(RVUNIFORM(0,1),d);        
rs = random(RV,N,1);
rs = cell2mat(rs);
fs = exp(-sum((rs-repmat(varargin{2},N,1)).^2*diag(varargin{3}.^2),2));         
%fs = prod(rs,2);
    case 5 %Borehole
RV = RANDVARS(RVNORMAL(0.1,0.0161812),RVNORMAL(0,1),...
    RVUNIFORM(63070,115600),RVUNIFORM(990,1110),RVUNIFORM(63.1,116),...
    RVUNIFORM(700,820),RVUNIFORM(1120,1680),RVUNIFORM(9855,12045));
rs = random(RV,N,1);
rs = cell2mat(rs);
fs = 2*pi*rs(:,3).*(rs(:,4)-rs(:,6))./(...
    log(exp(7.71+1.0056*rs(:,2))./rs(:,1)).*...
    (1+2*rs(:,7).*rs(:,3)./log(exp(7.71+1.0056*rs(:,2))./rs(:,1))./rs(:,1).^2./rs(:,8)+...
    rs(:,3)./rs(:,5))...
    ); 
RVxi = RANDVARS(RVNORMAL(),RVNORMAL(),RANDVARS(RVUNIFORM(-1,1),6));
rs = transfer(RV,RVxi,rs);
% fs=zeros(N,1);
% for k=1:N
% rw = xx(k,1);
% r  = exp(7.71+1.0056*xx(k,2));
% Tu = xx(k,3);
% Hu = xx(k,4);
% Tl = xx(k,5);
% Hl = xx(k,6);
% L  = xx(k,7);
% Kw = xx(k,8);
% 
% frac1 = 2 * pi * Tu * (Hu-Hl);
% 
% frac2a = 2*L*Tu / (log(r/rw)*rw^2*Kw);
% frac2b = Tu / Tl;
% frac2 = log(r/rw) * (1+frac2a+frac2b);
% 
% fs(k) = frac1 / frac2;
% end


end