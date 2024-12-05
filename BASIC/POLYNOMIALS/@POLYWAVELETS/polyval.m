function hval = polyval(h,K,x)
% hval = polyval(h,K,x)
%
% Evaluate a POLYWAVELETS h_n(x)
%
% x : vecteur (ou matrice) des points ou les polynomes sont evaluees
% K : liste des indices des polynomes a evaluer
% hval = matrice des valeurs . hval(i,j) = valeur de h_liste(j) en x(i)

param = get(h,'param');
n=param.n;
p=param.p;

switch min(size(x))
    case 0
        hval=sparse(0,length(K));
    case 1
        basis = @(x,r) double(randomeval(h.basis(r+1),x(:)))';
        mother = @(x,r) double(randomeval(h.mother(r+1),x(:)))';
        %generator = @(x) double((x(:)<1/2 & x(:)>=0)-(x(:)>=1/2 & x(:)<1));
        hval=zeros(numel(x),length(K));
        x=x(:);
        x=round(x*1e15)*1e-15;
        %pos = dec2bin(x*2^n);
        %xlocal = transfer(RVUNIFORM(0,1),RVUNIFORM(-1,1),x);
        level = zeros(getdim(h,p),1);
        %pos = zeros(getdim(h,p),1);
        level(1:(p+1))=-1;
        for k=0:n
            level((p+1)*(2^(k))+1:(p+1)*(2^(k+1)))=k;
        end
        Klevel = level(K+1);
        Kp = mod(K,p+1);
        repK = (Klevel==-1);
        hval(:,repK)=basis(x,Kp(repK));
        for k=0:n
            repK = find(Klevel==k);
            if ~isempty(repK)
                r = floor((K(repK)-2^k*(p+1))/(p+1));
                for m=0:2^k-1
                    repx = find( (x>=2^-k*(m) & x<2^-k*(m+1)) | (m==2^k-1 & x==1) );
                    hval(repx,repK(r==m))=2^(k/2)*mother(2^k*x(repx)-m,Kp(repK(r==m)));
                end
            end
            % if k==-1
            % elseif k==0 && K==1
            % hval(:,K==K)=2^(k/2)*generator(2^k*x);
            % else
            % l = K-(2^(k-1)-1)*(p+1)-2;
            % hval(:,K==K)=2^(k/2)*generator(2^k*x-l);
            % end
        end
        %hval = sparse(hval) ;
    otherwise
        hval = polyval(h,K,x(:));
        hval = reshape(full(hval),[size(x,1),size(x,2),length(K)]);
        
end