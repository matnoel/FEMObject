function c = expect(a,b)
% function c = expect(a,b)

switch nargin
    case 1
       c=mean(a); 
    case 2
       
       na = ndims(a.MYDOUBLE);
       nb = ndims(b.MYDOUBLE);
       a.MYDOUBLE = permute(a.MYDOUBLE,[a.pcdim,setdiff(1:na,a.pcdim)]);
       b.MYDOUBLE = permute(b.MYDOUBLE,[b.pcdim,setdiff(1:nb,b.pcdim)]);
       sa=size(a.MYDOUBLE);
       sb=size(b.MYDOUBLE);
       c = a.MYDOUBLE(:,:)' * b.MYDOUBLE(:,:) ;
       c = double(reshape(c,[sa(2:end),sb(2:end)]));
        
end
    