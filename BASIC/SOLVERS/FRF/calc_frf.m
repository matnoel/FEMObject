function [urw,uiw,umod,uphi]=calc_frf(scanw,K,C,M,f,ddloutput)
% function [urw,uiw,umod,uphi]=calc_frf(scanw,K,C,M,f,ddloutput)
% C : function handle: appel de C(w)
% f : function handle: appel de f(w)

N = size(K,1);
if isempty(C) 
    amort=0;
    if ~isa(f,'function_handle')
        b = @(w) f;
    else
        b = @(w) f(w);  
    end
    A = @(w) K-w^2*M;
else
    amort=1;
    if ~isa(f,'function_handle')
        b = @(w) [f;zeros(N,1)];
    else
        b = @(w) [f(w);zeros(N,1)];   
    end
    if ~isa(C,'function_handle')
        D = @(w) C;
    else
        D = @(w) C(w);
    end
    A = @(w) [K-w^2*M,-w*D(w);w*D(w),K-w^2*M] ;
end



urw = zeros(length(ddloutput),length(scanw));
uiw = zeros(length(ddloutput),length(scanw));

fprintf('Calcul de FRF : ')
for i=1:length(scanw)
    pourcentage(i,length(scanw),10)
    w=scanw(i);  
    u=solve(A(w),b(w));  
    ur = u(1:size(K,1));
    urw(:,i)=ur(ddloutput);    
    if amort
        ui = u(size(K,1)+1:end);    
        uiw(:,i)=ui(ddloutput);  
    end    
end

umod = sqrt(urw.^2+uiw.^2);
uphi = angle(urw+j*uiw);
